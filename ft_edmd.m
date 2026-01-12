function dataout = ft_edmd(cfg, datain)
% FT_EDMD  eDMD analysis integrated as a FieldTrip-style function

% Usage:
%   dataout = ft_edmd(cfg, datain)

% Inputs:
%   cfg - configuration struct (see defaults below)
%   datain - FieldTrip raw data struct (data.trial: 1 x Ntrials cell array)

% Outputs:
%   dataout - struct containing:
%     .cfg            - config used
%     .peakfeatures   - cell array (1 x Ntrials) of [Nsegments x Nbins] peak freqs
%     .sumfeatures    - cell array (1 x Ntrials) of [Nsegments x Nbins] power sums
%     .modefreqs      - cell array with per-trial/segment mode freqs (optional)
%     .modepowers     - cell array with per-trial/segment mode powers (optional)

%For loading ft_edmd first set up the cfg file: Set the defaults somewhere to look like this:
% cfg = [];
% cfg.D = 900;
% cfg.gamma = 200;
% cfg.nstacks = 5;
% cfg.seed = 1;
% cfg.trunk = 1e-4;
% cfg.freqEdges = [0 4 8 12 15 30 100];
% data_edmd = ft_edmd(cfg, data_eeg);  % where data is a FieldTrip raw data structure


% ----------------------------------
% FieldTrip Preamble  
% --------------------------
ft_revision = '$Id$';                                    %A Git placeholder
ft_nargin   = nargin;                                    %# of input arguments the user introduced
ft_nargout  = nargout;                                    %# of requested output arguments


% Ensure FieldTrip defaults are available
if exist('ft_defaults','file') == 2
    ft_defaults
else
    error('FieldTrip not found on path. Please don;t forget to run ft_defaults before calling ft_edmd.');
end

% preamble scripts (they run in caller workspace)
ft_preamble init                                          %runs ft_preamble_init.m inside the caller workspace.
ft_preamble debug                                         %debugging suppport
ft_preamble loadvar datain                                %override function inputs if the user requests loading from disk
ft_preamble provenance datain                             %attach provenance information about the input data

% ----------------------
% cfg defaults
% --------------------------

cfg.gamma    = ft_getopt(cfg, 'gamma', 200);
cfg.D        = ft_getopt(cfg, 'D', 900);
cfg.nstacks  = ft_getopt(cfg, 'nstacks', 5);
cfg.MA       = ft_getopt(cfg, 'MA', 0);
cfg.seed     = ft_getopt(cfg, 'seed', 1);
cfg.trunk    = ft_getopt(cfg, 'trunk', 1e-4);
cfg.freqEdges = ft_getopt(cfg, 'freqEdges', [0 4 8 12 15 30 100]);
cfg.verbose  = ft_getopt(cfg, 'verbose', true);
%cfg.channel = ft_getopt(cfg, 'channel', 'Fz', 'Cz', 'Pz'); %This is not properly formatted

% validate inputs
if ~isfield(datain, 'trial') || ~iscell(datain.trial)
    ft_error('datain must be a FieldTrip raw-like struct with datain.trial cell array.');
end

% binning parameters
freqEdges = cfg.freqEdges;
numBins = length(freqEdges)-1;

% set random seed for reproducibility
rng(cfg.seed);

% Preallocate outputs
nTrials = numel(datain.trial);
peakfeatures = cell(1, nTrials);
sumfeatures  = cell(1, nTrials);
modefreqs    = cell(1, nTrials);
modepowers   = cell(1, nTrials);

% ----------------------
% Main eDMD loop: (process each trial independently)
% ----------------------
for tr = 1:nTrials
    X_zero = datain.trial{tr};   % channels x time (assumed, as in any FT function)
    if isempty(X_zero)
        if cfg.verbose, ft_warning('trial %d is empty -> skipping', tr); end
        continue
    end

    % X_zero dimensions
    [nChannels, nSamples] = size(X_zero);

    % optional moving average smoothing along time
    if cfg.MA > 0
        w = ones(1, cfg.MA)/cfg.MA;
        for ch = 1:nChannels
            X_zero(ch, :) = conv(X_zero(ch, :), w, 'same');
        end
    end

    % Build stacked data (delay-embedding) if requested
    if cfg.nstacks > 1
        nst = cfg.nstacks;
        % create stacked matrix: (nChannels*nstacks) x newSamples
        newSamples = nSamples - nst + 1;
        Xaug = zeros(nChannels*nst, newSamples);
        for s = 1:nst
            Xaug((s-1)*nChannels + (1:nChannels), :) = X_zero(:, s:(s+newSamples-1));
        end
        X1 = Xaug(:, 1:end-1);
        X2 = Xaug(:, 2:end);
    else
        X1 = X_zero(:, 1:end-1);
        X2 = X_zero(:, 2:end);
    end

    % Update dims
    N = size(X1,1);   % dimension of stacked state
    M = size(X1,2);   % number of snapshots

    % --------------------------------------
    % RFF: random Fourier features (cos/sin pairs) dictionary
    % ----------------------
    D = cfg.D;
    W = randn(D, N) / cfg.gamma;      % D x N
    b = 2*pi*rand(D,1);               % D x 1  (random phase)
    WX1 = W * X1 + b * ones(1, size(X1,2));   % D x M
    WX2 = W * X2 + b * ones(1, size(X2,2));
    Psi_mX = sqrt(1/D) * [cos(WX1); sin(WX1)]; % (2D) x M
    Psi_mY = sqrt(1/D) * [cos(WX2); sin(WX2)];
    % transpose to M x (2D) for convenience in least-squares below if needed
    % but we'll use them in (2D x M) form for G/A calc as in your code

    % ----------------------
    % Build B / A / G / K
    % ----------------------
    Phi = Psi_mX.';      % M x (2D)    (note: earlier you used Phi = Psi_mX)
    % Fit linear map from observables -> each original channel via least squares
    % Acoeff: N x (2D)
%    Acoeff = zeros(N, size(Phi,2));
%    for ch = 1:N
%        y = X1(ch, :)';    % M x 1
%        a = Phi \ y;       % least squares solution (M x (2D) \ Mx1 -> (2D)x1)
%        Acoeff(ch, :) = a';
%    end
%here we terminate the least squares approach as is redundant
    B = pinv(Phi) * X1.';   % (2D) x N    <-- alternative full-state matrix (your original B)

    % eDMD (G, A, K)
    G = (Psi_mX * Psi_mX.') / M;   % (2D x 2D)
    A_mat = (Psi_mX * Psi_mY.') / M; % (2D x 2D)
    K = pinv(G) * A_mat;            % Koopman approx (2D x 2D)

    % ----------------------
    % SVD and truncation (economy)
    % ----------------------
    try
        [U,S,V] = svd(K, 'econ');
    catch ME
        warning('SVD on K failed for trial %d: %s', tr, ME.message);
        U = []; S = []; V = [];
    end

    diagEntries = diag(S);
    r = find(diagEntries >= cfg.trunk, 1, 'last');
    if isempty(r)
        r = min(10, length(diagEntries)); % fallback minimal rank
    end
    U_t = U(:, 1:r);
    S_t = S(1:r, 1:r);
    V_t = V(:, 1:r);
    Kt = U_t*S_t*V_t';

    % ----------------------
    % Eigendecomposition of truncated K
    % ----------------------
    [Wleft, MU] = eig(Kt);  % Wleft contains right eigenvectors? keep naming consistent with your pipeline
    % Note: you might want left vs right eigenvectors depending on mode calc; adjust if needed
    % MU is diagonal matrix with eigenvalues

    % compute normalization of left/right eigenvectors (adapting your steps)
    XI = []; % placeholder for right eigenvectors if needed
    try
        % If XI (right eigenvectors) needed: get them from eig on Kt' or from V
        [XI_tmp, ~] = eig(Kt'); % right eigenvectors of Kt' are left eigenvectors of Kt
        XI = XI_tmp;
    catch
        XI = [];
    end

    % some of your normalization steps (safe-guarded)
    try
        W = Wleft; % align with your notation
        if ~isempty(XI)
            inerp = sum(conj(W).*XI);
            Wn = W ./ inerp;
            scaledIP = sum(conj(Wn).*XI);
            ip_angle = angle(scaledIP);
            an_cor = exp(1i * ip_angle);
            Wnn = Wn .* an_cor;
        else
            Wnn = W;
        end
    catch
        Wnn = Wleft;
    end

    % Koopman modes (using your formula V = (Wnn'*B).')
    try
        V_modes = (Wnn' * B).';   % N x r
        V_modes(isnan(V_modes)) = 0;
    catch
        V_modes = zeros(N, size(Wnn,2));
    end

    % ----------------------
    % State reconstruction (optional & approximate)
    % ----------------------
    % Form discrete-time eigenvalues and their temporal evolution
    mu_diag = diag(MU).';   % 1 x r
    % build mu^k matrix for k = 1..M (we will use powers for reconstruction)
    row_idx = (1:M).';
    mu_p = mu_diag .^ row_idx;   % M x r  (elementwise)
    % use initial observable Psi0 (first column of Psi_mX)
    Psi0 = Psi_mX(:,1).';  % 1 x (2D)
    % attempt phi diag using XI (if available)
    if ~isempty(XI)
        Phid = diag(Psi0 * XI);
        Xrr = real( (mu_p * Phid * V_modes')' );  % N x M
    else
        Xrr = nan(size(X1)); % skip if unavailable
    end

    % Reduce back to original number of channels if stacking was used
    orN = size(X_zero,1);
    if size(Xrr,1) >= orN
        Xr_rec = Xrr(1:orN, :);
    else
        Xr_rec = nan(orN, size(Xrr,2));
    end

    % ----------------------
    % Frequency and power calculation per mode
    % ----------------------
    % compute power per mode (norm squared)
    P_modes = vecnorm(V_modes).^2;
    % convert discrete eigenvalues MU to frequency estimates
    % typical mapping: freq = angle(lambda)/(2*pi*dt)  ; dt unknown — we assume dt=1 or scale externally
    % Use your empirical scaling constant if desired; provide raw freq (cycles per sample)
    dt = 1; % user may change or pass in sampling freq if available (datain.fsample)
    f_modes = abs(angle(diag(MU)) / (2*pi*dt));
    % threshold tiny freqs
    f_modes(f_modes < 1e-6) = 0;

    % Filter out zero-frequency modes (e.g., stationary)
    validIdx = f_modes > 0;
    fValid = f_modes(validIdx);
    PValid = P_modes(validIdx);

    % Sort & unique similar to your code
    [fSorted, I] = sort(fValid);
    PSorted = PValid(I);
    [~, uniqueIdx] = unique(fSorted, 'stable');
    fUnique = fSorted(uniqueIdx);
    PUnique = PSorted(uniqueIdx);

    % Binning: compute peak frequency and sum power per bin
    peakvec = computeBinnedPeaks(fUnique, PUnique, freqEdges);
    sumvec  = computePowerSum(fUnique, PUnique, freqEdges);

    % Save outputs per trial
    peakfeatures{tr} = peakvec;
    sumfeatures{tr}  = sumvec;
    modefreqs{tr}    = fUnique;
    modepowers{tr}   = PUnique;

    if cfg.verbose
        fprintf('ft_edmd: trial %d processed (modes: %d -> kept %d unique freqs)\n', tr, numel(f_modes), numel(fUnique));
    end
end

% ----------------------
% Wrap outputs into FT-style struct
% ----------------------
dataout = [];
dataout.cfg = cfg;
dataout.peakfeatures = peakfeatures;
dataout.sumfeatures = sumfeatures;
dataout.modefreqs = modefreqs;
dataout.modepowers = modepowers;

% ----------------------
% Postamble (FieldTrip style)
% ----------------------
ft_postamble debug
ft_postamble previous datain
ft_postamble provenance dataout
ft_postamble history dataout
ft_postamble savevar dataout

end

% ----------------------
% Helper functions (local)
% ----------------------
function innerbinpeaks = computeBinnedPeaks(f, P, freqEdges)
    nb = length(freqEdges)-1;
    innerbinpeaks = NaN(1, nb);
    for b = 1:nb
        idx = f >= freqEdges(b) & f < freqEdges(b+1);
        localpower = P(idx);
        localfreq  = f(idx);
        if isempty(localpower)
            innerbinpeaks(b) = NaN;
        else
            [~, Didx] = max(localpower);
            innerbinpeaks(b) = localfreq(Didx);
        end
    end
end

function powersum = computePowerSum(f, P, freqEdges)
    nb = length(freqEdges)-1;
    powersum = NaN(1, nb);
    for b = 1:nb
        idx = f >= freqEdges(b) & f < freqEdges(b+1);
        localpower = P(idx);
        if isempty(localpower)
            powersum(b) = NaN;
        else
            powersum(b) = sum(localpower);
        end
    end
end
