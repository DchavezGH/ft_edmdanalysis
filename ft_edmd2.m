function dataout = ft_edmd2(cfg, datain)
% FT_EDMD  performs eDMD decomposition on time series data

% Usage:
%   dataout = ft_edmd(cfg, datain)

%   Mandatory inputs:
%   cfg.reconstruct = 'yes' or 'no', compute state reconstruction Xr_rec (default = 'no')

%   Optional inputs: 
%   cfg         -  Configuration struct (see defaults below)
%   datain      -  FieldTrip raw data struct (data.trial: 1 x Ntrials cell array)
%   cfg.cut        = double, rank truncation threshold (default = 0.9999999999)
%   cfg.gamma      = double, RFF scaling parameter (default = 4)
%   cfg.D          = double, number of random Fourier features (default = 900)
%   cfg.nstacks    = double, Hankel stacking depth (default = 5)
%   cfg.MA         = double, moving average window (default = 0)
%   cfg.seed       = double, RNG seed (default = 1)
%   cfg.trunk      = double, SVD truncation tolerance (default = 1e-4)
%   cfg.freqEdges  = double vector, frequency bin edges (default = [0 4 8 12 15 30 100])
%   cfg.verbose    = logical, print progress (default = true)

% The configuration can also contain any parameters needed for the specific
% aanalysis:

%   cfg.channel
%   cfg.trials
%   cfg.latency 
%   etc...

% Outputs:
%   dataout - struct containing:
%     .cfg            - config used
%     .peakfeatures   - cell array (1 x Ntrials) of [Nsegments x Nbins] peak freqs
%     .sumfeatures    - cell array (1 x Ntrials) of [Nsegments x Nbins] power sums
%     .modefreqs      - cell array with per-trial/segment mode freqs (optional)
%     .modepowers     - cell array with per-trial/segment mode powers (optional)


%% -----------------------------
%% Validate input data
%% -----------------------------
datain = ft_checkdata(datain, ...
    'datatype', {'raw','raw+comp'}, ...
    'hassampleinfo', 'yes');

cfg = ft_checkconfig(cfg, 'required', {});


% -----------------------------
% Set configuration defaults
% -----------------------------
cfg = ft_checkconfig(cfg, 'deprecated', {}); %Related to backwards compatibility
cfg = ft_checkconfig(cfg, 'renamed', {});    %Related to backwards compatibility


%cfg = ft_checkopt(cfg, 'D', 'double');
cfg = ft_checkopt(cfg, 'reconstruct', 'char', {'yes','no'});

cfg.reconstruct = ft_getopt(cfg, 'reconstruct', 'no');
cfg.cut         = ft_getopt(cfg, 'cut', 0.9999999999);
cfg.gamma       = ft_getopt(cfg, 'gamma', 4);
cfg.D           = ft_getopt(cfg, 'D', 900);
cfg.nstacks     = ft_getopt(cfg, 'nstacks', 5);
cfg.MA          = ft_getopt(cfg, 'MA', 0);
cfg.seed        = ft_getopt(cfg, 'seed', 1);
cfg.trunk       = ft_getopt(cfg, 'trunk', 1e-4);
cfg.freqEdges   = ft_getopt(cfg, 'freqEdges', [0 4 8 12 15 30 100]);
cfg.verbose     = ft_getopt(cfg, 'verbose', true);

%% -----------------------------
%% Channel / trial / latency selection
%% -----------------------------
% Only use channel selection if field exists
if isfield(cfg,'channel')
    % Only keep channels that actually exist
    cfg.channel = intersect(cfg.channel, datain.label);
end

datain = ft_selectdata(cfg, datain);

% Validate again after selection
if isempty(datain.trial)
    ft_error('No trials left after selection. Check your cfg.channel, cfg.trials, and cfg.latency.');
end

%% -----------------------------
%% Binning parameters
%% -----------------------------
freqEdges = cfg.freqEdges;
numBins   = length(freqEdges)-1;

%% -----------------------------
%% Seed RNG for reproducibility
%% -----------------------------
%rng(cfg.seed);
rng(cfg.seed, 'twister');

%% -----------------------------
%% Preallocate outputs
%% -----------------------------
nTrials      = numel(datain.trial);
peakfeatures = cell(1, nTrials);
sumfeatures  = cell(1, nTrials);
modefreqs    = cell(1, nTrials);
modepowers   = cell(1, nTrials);
%xrecon       = cell(1, nTrials);
if strcmp(cfg.reconstruct, 'yes')
    xrecon = cell(1, nTrials);
end



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

    % Build stacked data (Henkel-matrix)
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

    % ------------------------------------
    % RFF: random Fourier features (cos/sin pairs) dictionary
    % ------------------------------------
    D = cfg.D;
    W = randn(D, N) / cfg.gamma;      % D x N
    b = 2*pi*rand(D,1);               % D x 1  (random phase)
    WX1 = W * X1 + b * ones(1, size(X1,2));   % D x M
    WX2 = W * X2 + b * ones(1, size(X2,2));
    Psi_mX = (sqrt(1/D) * [cos(WX1); sin(WX1)])'; % (2D) x M
    Psi_mY = (sqrt(1/D) * [cos(WX2); sin(WX2)])';

    % ----------------------
    % Build compact K operator: Phi / B / U / S / V / Kt2 / 
    % ----------------------
      
    [U,S,V] = svd(Psi_mX / sqrt(M), 'econ');  
    singvals = diag(S);
    r = find(cumsum(singvals) >= cfg.cut * sum(singvals), 1, 'first');
    U = U(:,1:r);
    S = S(1:r,1:r);
    V = V(:,1:r);
    Kt2 = (S \ (U' * (Psi_mY/sqrt(M)) * V));    
    B = V*pinv(S)*U' * X1';   %using SVD of Psi_mX (new, fast full-state matrix:)
    
    % ----------------------
    % Eigendecomposition of K: Wnn / MU / XI
    % ----------------------
    [XI,MU,W] = eig(Kt2);       %Generate eigenvalues and left-right eigenvctors
    inerp = sum(conj(W).*XI);                  %Calculate the w_k'xi_k products
    %inerp(inerp< 0.00001) = 0;   %Remove neglible inner p. (0.00001 worked ok) 
    Wn = W./inerp;    %define the scaled w_n left eigenvectors in the u. circle
    scaledIP = sum(conj(Wn).*XI);              %Calculate the w_n'xi_k products
    ip_angle = angle(scaledIP);     %phase of w_n'xi_k products in the c. plane
    an_cor = exp(1i * ip_angle);                      %angle correction factor
    Wnn = Wn.*an_cor;                                 %normalized w eigevectors
    %ipfinal = sum(conj(Wnn).*XI); %this should be vector of ones if all went ok

    % ----------------------
    % Koompan modes, and state reconstruction: V_modes  / MU / XI
    % ----------------------

    %Calculate the Koopman modes v_i
    V_modes = (Wnn'*V'*B).';    %We must project l-eigenvectors on V' to recover time dimensions
    V_modes(isnan(V_modes)) = 0; %remove all undefined entries.
    
if strcmp(cfg.reconstruct, 'yes')
    mu = repmat(diag(MU)', M, 1);
    row_indices = (1:M)';  % Create a column vector of row indices
    mu_p = mu .^ row_indices;  % Element-wise exponentiation
    Psi0 = Psi_mX(1,:);
    %State reconstruction: build the diagonal Phi matrix
    Phid = diag(Psi0*V*XI);  %We must project r-eigenvectors on V to recover time dimensions
    Xrr = real((mu_p*Phid*V_modes')'); %State recontruction: 
    orN = nChannels; %Original number of channels nChannels
    Xr = Xrr(1:orN,:); %Remove copies from the stacking 
    Xr_col = Xr(:);
    Xr_col_n = normalize(Xr_col, "range");
    Xr_rec = reshape(Xr_col_n, orN, M);
end
    % ----------------------
    % Frequency and power calculation per mode
    % ----------------------
    % compute power per mode (norm squared)
    P_modes = vecnorm(V_modes).^2;
    % convert discrete eigenvalues MU to frequency estimates
    %dt = 1; % user may change or pass in sampling freq if available (datain.fsample)
    dt = 1 / datain.fsample;
    f_modes = abs(angle(diag(MU)) / (2*pi*dt));
    f_modes(f_modes < 1e-6) = 0;     % threshold tiny freqs

    % Filter out zero-frequency modes 
    validIdx = f_modes > 0;
    fValid = f_modes(validIdx);
    PValid = P_modes(validIdx);

    % Sort & unique
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
  if strcmp(cfg.reconstruct, 'yes')
    xrecon{tr}       = Xr_rec;
  end
    ft_info('trial %d processed ...', tr);
    %if cfg.verbose
        %ft_info('ft_edmd: trial %d processed (modes: %d -> kept %d unique freqs)\n', tr, numel(f_modes), numel(fUnique));
    %end
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
if strcmp(cfg.reconstruct, 'yes')
dataout.xrecon = xrecon;
end
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
% Helper functions
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
