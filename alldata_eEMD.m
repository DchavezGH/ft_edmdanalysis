%This code takes all subjects, all trials, and performs eDMD on them

freqEdges =[0 4 8 12 15 30 100];  %neuro-informed bin edges
numBins = length(freqEdges) - 1;
sumfeats = cell(numSubjects, runsPerSubject); %Output cells 
peakfeats = cell(numSubjects, runsPerSubject); %Output cells

for subj = 1:numSubjects

    for run = 1:runsPerSubject
        if ~isempty(fmidData{subj, run})
            disp(['Processing cell ' num2str(subj)])
            pSegmentstotal = fmidData{subj, run};  %Inner name for the data
            numSegments = numel(pSegmentstotal); %pSegments cardinality
            sumfeatures = zeros(numSegments, numBins);
            peakfeatures = zeros(numSegments, numBins);
            
            for k = 1:numSegments 
       %    for k = [3 ]
    

    
            X_zero = pSegmentstotal{k};        % matrix for this segment
            % run eDMD:
            %eDMD parameters
            %gamma = 1.3;  %Standard value, used for all previous analysis 
            gamma = 200; %Quite ok

            nstacks = 5;                                          %Stacking parameter
            N = size(X_zero,1);                 %ARTIFICIAL DATASET: number of channels in the data
            M = size(X_zero,1);              %ARTIFICIAL DATASET: number of snapshots in the data
            t = linspace(0, 1, M); %This is needed further down to build B
            MA = 3;                      %ARTIFICIAL DATASET: Moving average parameter
            D =900; %600;              % Half the number of total features (final dim = 2D)
            X1 = X_zero(:, 1:end-1);
            X2 = X_zero(:, 2:end);


            %Here we will perform stacking: 
            if nstacks > 1
                Xaug = [];
                for st = 1:nstacks
                    Xaug = [Xaug; X_zero(:, st:end-nstacks+st)]; 
                end
                
                X1 = Xaug(:, 1:end-1);
                X2 = Xaug(:, 2:end);
            else
                X1 = X_zero(:, 1:end-1);
                X2 = X_zero(:, 2:end);
            end



            %Stacked parameters. No need to define them in advance, but they need to be
            %redefined to reflect the new dimentions of the dataset:
            N = size(X1,1);  %number of channels in the stack
            M = size(X1,2); %number of snapshots in the stack 
            t = linspace(0, 1, M);

            %rng(18347);               % Seed for reproducibility
           rng(1);               

            %RFF SETUP 
            N = size(X1, 1);                      % Input dimension
            W = randn(D, N) / gamma;             % Random frequency vectors: D x N with variance (1/gamma)^2
%            b = 2 * pi * rand(D, 1);             % Random phase shifts: D x 1
            b = 0 * pi * rand(D, 1);             % Random phase shifts: D x 1
   
            % RFF transformation with cos/sin pairs
            WX1 = W * X1 + b * ones(1, size(X1,2));  % Now D×M
            WX2 = W * X2 + b * ones(1, size(X2,2));
            
            Psi_mX = sqrt(1 / D) * [cos(WX1); sin(WX1)];  % (2D) x M
            Psi_mY = sqrt(1 / D) * [cos(WX2); sin(WX2)];
            
            Psi_mX = Psi_mX'; %These are the obserbable matrices from the paper 
            Psi_mY = Psi_mY'; %These are the obserbable matrices from the paper 

            %GENERATE THE B MATRIX: weights that are used to build the koopman modes

            % Build design matrix B
            Phi = Psi_mX;
            %Phi should actually contain the elements of Psi_mY, so the least squares
            %problem delivers a useful coeff. matrix: 

            % Preallocate coefficients
            A = zeros(N, size(Phi,2));   
            % Fit each channel's time-series data using least squares
            for m = 1:N
                y = X1(m, :)';                 % Column vector for time series of channel m
                a = Phi \ y;                  % Solve least squares: Phi * a ≈ y
                A(m, :) = a';                 % Store in row
            end
            lqrec = Phi * A';
            B=pinv(Phi)*X1'; %Full state matrix, note we don't need identity in the dictionary
            Run wthe next couple lines to verify that B fulfils eq16 in Williams:
            % g=B.'*Psi_mX.';   %the full state observable g(x)
            % CC=g-X1;  %error between g(x) and original data
            G = (Psi_mX'*Psi_mX)/M;
            A = (Psi_mX'*Psi_mY)/M;
            K = pinv(G)*A;    %The eDMD approximation to the Koopman operator
            %clear G A

            %We introduce SVD to truncate K
            [u,sig,v] = svd(K, 'econ'); % Performs svd on the data

            %automatic Truncation deletes eigs<0.0001
            trunk = 0.0001; %eigenvalues smaller than 0.0001
            diagEntries = diag(sig);  % extract the diagonal entries as a vector
            r = find(diagEntries >= trunk, 1, 'last') + 1; % First index where all subsequent are < trunk         
            %SVD-truncated K operator
            u_t = u(:,1:r); 
            sig_t = sig(1:r,1:r); 
            v_t = v(:,1:r); % r-truncate U, Sigma & V
            Kt = u_t*sig_t*v_t' ;   %Generates the truncated Kt

            %Calculate and normalize eigenvectors: 
            [XI,MU,W] = eig(Kt,"nobalance");       %Generate eigenvalues and left-right eigenvctors
            %[XIb,MUb,Wb] = eig(Kt); %balanced eigenvalues
            inerp = sum(conj(W).*XI);                  %Calculate the w_k'xi_k products
            Wn = W./inerp;    %define the scaled w_n left eigenvectors in the u. circle
            scaledIP = sum(conj(Wn).*XI);              %Calculate the w_n'xi_k products
            ip_angle = angle(scaledIP);     %phase of w_n'xi_k products in the c. plane
            an_cor = exp(1i * ip_angle);                      %angle correction factor
            Wnn = Wn.*an_cor;                                 %normalized w eigevectors
            ipfinal = sum(conj(Wnn).*XI); %this should be vector of ones if all went ok

            %Calculate the Koopman modes v_i
            V = (Wnn'*B).';    %eDMDSVD: We use conjugate transpose now...
            V(isnan(V)) = 0; %remove all undefined entries.

            %State reconstruction with discrete-time eigenvalues:
            mu = repmat(diag(MU)', M, 1);
            row_indices = (1:M)';  % Create a column vector of row indices
            mu_p = mu .^ row_indices;  % Element-wise exponentiation
            Psi0 = Psi_mX(1,:);
            Phid = diag(Psi0*XI);  %State reconstruction: build the diagonal Phi matrix
            Xrr = real((mu_p*Phid*V')'); %State recontruction: 
            orN = size(X_zero,1); %Original number of channels
            Xr = Xrr(1:orN,:); %Remove copies from the stacking 
            Xr_col = Xr(:);
            Xr_col_n = normalize(Xr_col, "range");
            Xr_rec = reshape(Xr_col_n, orN, M);


            %Freq and Power 
            P = vecnorm(V).^2; % power per mode
            %f = 9.0746.*diag(abs(imag(MU) ./ (2*pi*0.004))); %added a correction factor 
            f = 1.444.*diag(abs(imag(MU) ./ (2*pi*0.004))); %new correction factor (f/2*pi)
            f(f < 0.01) = 0; % threshold to remove small freqs
            
            % Remove zero or duplicate frequencies
            validIdx = f > 0;           % find zero-frequency modes
            fValid = f(validIdx);
            PValid = P(validIdx);   % drop zero-frequency modes
            [fSorted, I] = sort(fValid);
            PSorted = PValid(I);
            [~, uniqueIdx] = unique(fSorted, 'stable'); 
            fUnique = fSorted(uniqueIdx);
            PUnique = PSorted(uniqueIdx);
            
            % Now safe to run findpeaks
            peakfeatures(k, :) = computeBinnedPeaks(fUnique, PUnique, freqEdges);
            sumfeatures (k, :) = computePowerSum(fUnique, PUnique, freqEdges);
            %from inner variable to output cell:
  

            %GRAPHS GRAPHS GRAPHS GRAPHS GRAPHS GRAPHS GRAPHS GRAPHS GRAPHS
        % % %Graphs: SVD singular values (and truncated sv)
        % figure;
        % subplot(1,2,1)
        % plot(sig)
        % title('SVD SIGMA');
        % subplot(1,2,2)
        % plot(sig_t)
        % title('truncated SVD sigma');


        % % %Graphs power and frequency per mode
        % figure;
        % plot(fUnique,PUnique)
        % title('freq vs power')

        % %Graphs: Artificial dataset and eDMD reconstruction
        % figure;
        % subplot(2,1,1)
        % surf(X_zero, 'EdgeColor','none')
        % title('artificial dataset', gamma)
        % subplot(2,1,2)
        % surf(Xr_rec, 'EdgeColor','none')
        % title('eDMD reconstruction')

%         % %PRINTABLE VERSION START: 
% figure;
% lineColor = [0.5 0.5 0.5]; % grey lines
%   subplot(2,1,1)  %----SUBPLOT 1------
% 
%  surface
% surf(X_zero,'FaceAlpha',1);
% %title('artificial dataset', k)
% shading interp;
% %campos([587.72, -25.9, 7.5]);
% hold on
% for r = [1 2 3 4 5 6 7]
%     plot3(1:size(X_zero,2), ...         % x-axis
%           r*ones(1,size(X_zero,2)), ... % y-axis (row index)
%           X_zero(r,:), ...              % z values
%           'Color', lineColor, 'LineWidth', 0.35);
% end
% hold off
% % Fonts
% ax = gca;
% % ax.CameraPosition = [587, -25.9, 7.5];
% % ax.CameraTarget = [35, 35, 0];
% %ax.CameraViewAngle = 6;          
% ax.YTick = [36];
% ax.ZTick = [-100 100];
% ax.FontName = 'Helvetica';
% ax.FontSize = 10;
% 
% %title('a', 'FontName', 'Helvetica', 'FontSize', 12);
% text(ax, 0.9, 1, 'a)', 'Units', 'normalized', ...
%      'FontName', 'Helvetica', 'FontSize', 11, 'FontWeight', 'normal');
% 
% colormap(parula);   % or viridis/cbrewer if installed
% 
% subplot(2,1,2)   %----SUBPLOT 2------
%  surface
% surf(Xr_rec,'FaceAlpha',1);
% caxis([0.1*max(Xr_rec(:)) 1]);
% shading interp;
% 
% hold on
% for r = [1 2 3 4 5 6 7]
%     plot3(1:size(Xr_rec,2), ...
%           r*ones(1,size(Xr_rec,2)), ...
%           Xr_rec(r,:), ...
%           'Color', lineColor, 'LineWidth', 0.35);
% end
% hold off
% % Fonts
% ax = gca;
% % ax.CameraPosition = [587, -25.9, 7.5];
% % ax.CameraTarget = [35, 35, 0];
% % ax.CameraViewAngle = 6;   
% ax.YTick = [36];
% ax.ZTick = [0 1];
% ax.FontName = 'Helvetica';
% ax.FontSize = 10;
% 
% %title('b', 'FontName', 'Helvetica', 'FontSize', 10);
% text(ax, 0.9, 1, 'b)', 'Units', 'normalized', ...
%      'FontName', 'Helvetica', 'FontSize', 11, 'FontWeight', 'normal');
% colormap(parula);   % or viridis/cbrewer if installed
% 
% 
% 
% % Line width for mesh edges (if any)
% set(gca, 'LineWidth', 0.75);
% 
% % Set figure size: 76 mm wide
% %width = 76/10; height = 60/10;
% width = 76/10; height = 60/10;
% set(gcf, 'Units', 'centimeters', 'Position', [5 5 width height]);
% set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [width height]);
% 
% % Export as PDF
% 
%  filename = sprintf('surface_example_subj%d_run%d.pdf', subj, k);
%  print(gcf, filename, '-dpdf', '-painters');
% 
% 
% %%PRINTABLE VERSION END

% %PRINTABLE VERSION OF THE FREQ PLOT: 
% figure;
% 
% % Plot with blue line, 0.75 pt thickness
% plot(fUnique, PUnique, 'b-', 'LineWidth', 0.75);
% 
% % Title and axes with specified fonts and sizes
% %title('freq vs power', 'FontName', 'Helvetica', 'FontSize', 12);
% xlabel('Frequency', 'FontName', 'Helvetica', 'FontSize', 10);
% ylabel('Power', 'FontName', 'Helvetica', 'FontSize', 10);
% 
% % Set axis tick font
% ax = gca;
% ax.FontName = 'Helvetica';
% ax.FontSize = 10;
% ax.XTick = [4 8 12 20 30 45];
% % Set figure size: 76 mm wide, keep aspect ratio (height can be adjusted)
% width = 76/10;   % mm → cm (MATLAB expects cm)
% height = 40/10;  % choose proportion, e.g. 60 mm here
% set(gcf, 'Units', 'centimeters', 'Position', [5 5 width height]);
% set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [width height]);
% 
% % Export as vector PDF
% print(gcf, 'fvp33_noise.pdf', '-dpdf', '-painters');
% %PRINTABLE VERSION END









end
            sumfeats{subj, run} = peakfeatures; 
            peakfeats{subj, run} = sumfeatures;


            

        end
    end
end


















%FUNCTIONS

function binnedPower = computeBinnedPower(f, P, freqEdges) %chat's idea to bin power:
    % Initialize binned power vector
    binnedPower = zeros(1, length(freqEdges)-1);
    
    % For each bin, sum power of frequencies inside bin range
    for b = 1:length(freqEdges)-1
        idx = f >= freqEdges(b) & f < freqEdges(b+1);
        binnedPower(b) = sum(P(idx));  %This is the part we needa modify
        %[pks, locs] = findpeaks(P, f, NPeaks=1, SortStr="descend");
    end
end

function innerbinpeaks = computeBinnedPeaks(f, P, freqEdges) %david's idea to bin power:
    % Initialize binned power vector
    innerbinpeaks = zeros(1, length(freqEdges)-1);
    
    % For each bin, find max power value of frequencies inside bin range
    for b = 1:length(freqEdges)-1
                % idx = f >= freqEdges(b) & f < freqEdges(b+1); %logic array 0 if condition is not met, else 1
                % %P(idx)
                % [pks, locs] = findpeaks(P(idx), f(idx), NPeaks=1);
                % innerbinpeaks(b)=locs
                
                idx = f >= freqEdges(b) & f < freqEdges(b+1); %logic array 0 if condition is not met, else 1
                localpower = P(idx);
                localfreq = f(idx);
                %[pks, Didx] = max(localpower);
                %innerbinpeaks(b)= localfreq(Didx);
                if isempty(localpower)
                   innerbinpeaks(b) = NaN; % No value for this bin
                else
                   [~, Didx] = max(localpower);
                   innerbinpeaks(b) = localfreq(Didx);
                end
                
    end
end

function powersum = computePowerSum(f, P, freqEdges) %david's idea to bin power:
    % Initialize binned power vector
    powersum = zeros(1, length(freqEdges)-1);

    % For each bin, sum power of frequencies inside bin range
    for b = 1:length(freqEdges)-1
                % idx = f >= freqEdges(b) & f < freqEdges(b+1); %logic array 0 if condition is not met, else 1
                % %P(idx)
                % [pks, locs] = findpeaks(P(idx), f(idx), NPeaks=1);
                % innerbinpeaks(b)=locs

                idx = f >= freqEdges(b) & f < freqEdges(b+1); %logic array 0 if condition is not met, else 1
                localpower = P(idx);
                localfreq = f(idx);
                %[pks, Didx] = max(localpower);
                %innerbinpeaks(b)= localfreq(Didx);
                if isempty(localpower)
            powersum(b) = NaN; % No value for this bin
        else
            powersum(b) = sum(localpower); 
        end

    end
end
