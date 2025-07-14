%This code generates an artificial dataset (with  two oscillatory signals),
% performs eDMD on it via a random Fourier features dictionary,and generates
% power VS frequency graphs. Then, it cycles through values of gamma (from 
% the RFF) dictionary so an "optimal value" can be determined


%for siggma = 1.0:-0.2:0.2
%0.6 and 1 kinda work
%for siggma = 0.7:-0.02:0.5
%0.54 0.60 0.62 0.64 0.66 
%for siggma = 0.1:0.01:0.2
%0.11 works
%for siggma = 0.10:0.001:0.12
%0.11 works haha
%for siggma = 0.0010:0.00001:0.0012
%0.00113 kinda works
%for siggma = [0.00113, 0.11, 0.6, 1]
siggma = [0.0113];

nstacks = 5;                                          %Stacking parameter
r = 152;                                              %SVD: rank truncation 
N = 2;                 %ARTIFICIAL DATASET: number of channels in the data
M = 160;              %ARTIFICIAL DATASET: number of snapshots in the data
%t = dt.*(1:M);       %ARTIFICIAL DATASET: Time vector, temporarily ignored
t = linspace(0, 1, M); %This is needed further down to build B
omega = 17.5;                      % ARTIFICIAL DATASET: Frequency multiplier
tau = 1.2;                            % ARTIFICIAL DATASET: delay parameter
MA = 3;                      %ARTIFICIAL DATASET: Moving average parameter
D = 600;              % Half the number of total features (final dim = 2D)


%ARTIFICIAL DATASET:
Nn= N;     %Safekeeping original channel dim
%Xs = -4.*exp(-(12.*t-7).^2).*sin(100.*t); %This works
%Xss = 0.9.*sin(60.*t);    %This works
Xs = -4.*exp(-(12.*t-7).^2).*sin(67.*t); 
Xss = 0.2.*sin(84.*t);    
Xi = repmat(Xs, N/2, 1);
Xii = repmat(Xss, N/2, 1);
X_zero = [Xi;Xii];
%X_zero = randn(size(X_zero ));
X_0_rec = X_zero; %take it off when retaking normalization

X1 = X_zero(:, 1:end-1);
X2 = X_zero(:, 2:end);


%Here we will perform stacking: 
%nstacks = 1;  %Stacking, at the parameters section
if nstacks > 1
    Xaug = [];
    for st = 1:nstacks
        Xaug = [Xaug; X_0_rec(:, st:end-nstacks+st)]; 
    end
    
    X1 = Xaug(:, 1:end-1);
    X2 = Xaug(:, 2:end);
else
    X1 = X_0_rec(:, 1:end-1);
    X2 = X_0_rec(:, 2:end);
end



%Stacked parameters. No need to define them in advance, but they need to be
%redefined to reflect the new dimentions of the dataset:
N = size(X1,1);  %number of channels in the stack
M = size(X1,2); %number of snapshots in the stack 
t = linspace(0, 1, M);

%Dictionary
% X1e = [ X1.^2 ; X1.^3 ; X1.^4 ; X1.^5 ; X1.^6 ; X1.^7 ; X1.^8 ; X1.^9];
% X2e = [ X2.^2 ; X2.^3 ; X2.^4 ; X2.^5 ; X2.^6 ; X2.^7 ; X2.^8 ; X2.^9]; 
% 
% Psi_mX = X1e.'; %Proper variable naming
% Psi_mY = X2e.';



%D = 600;              % Half the number of total features (final dim = 2D)
sigma = siggma;          % Bandwidth for the RBF kernel
rng(1);               % Seed for reproducibility

% --- RFF MAPPING SETUP ---
N = size(X1, 1);                      % Input dimension
W = randn(D, N) / sigma;             % Random frequency vectors: D x N
b = 2 * pi * rand(D, 1);             % Random phase shifts: D x 1

% --- RFF TRANSFORMATION WITH COS/SIN PAIRS ---
WX1 = W * X1 + b * ones(1, size(X1,2));  % Now D×M
WX2 = W * X2 + b * ones(1, size(X2,2));

Psi_mX = sqrt(1 / D) * [cos(WX1); sin(WX1)];  % (2D) x M
Psi_mY = sqrt(1 / D) * [cos(WX2); sin(WX2)];

Psi_mX = Psi_mX'; 
Psi_mY = Psi_mY';

%GENERATE THE B MATRIX: weights that are used to build the koopman modes

%degree = 9; %Dictionary degree DELETE
%t_stack = dt.*(1:M);   %t vector for the dims of the stack

% Build design matrix 
%Phi = Psi_mX';
Phi = Psi_mX;
%Phi should actually contain the elements of Psi_mY, so the least squares
%problem delivers a useful coeff. matrix: 

%January mod: take away A
% Preallocate coefficients
% %A = zeros(N, (degree-1)*(N)); %Rows: channels, Columns: polynomial terms DELETE
% A = zeros(N, size(Phi,2));   
% % Fit each channel's time-series data using least squares
% for m = 1:N
%     y = X1(m, :)';                 % Column vector for time series of channel m
%     a = Phi \ y;                  % Solve least squares: Phi * a ≈ y
%     A(m, :) = a';                 % Store in row
% end
% 
% 
% lqrec = Phi * A';

B=pinv(Phi)*X1';
% g=B.'*Psi_mX.';


%Run the next couple lines to verify
%that B fulfils eq16 in Williams:
g=B.'*Psi_mX.';   %the full state observable g(x)
%g=B.'*Psi_mX;   %the full state observable g(x)
CC=g-X1;  %error between g(x) and original data

%January mod: take away GAK:
G = (Psi_mX'*Psi_mX)/M;  
A = (Psi_mX'*Psi_mY)/M;
K = pinv(G)*A;    %The eDMD approximation to the Koopman operator
% clear G A
%We introduce SVD to truncate K
[u,sig,v] = svd(K, 'econ'); % Performs svd on the data
%r=5; %rank truncation parameter is at the parameters section
u_t = u(:,1:r); 
sig_t = sig(1:r,1:r); 
v_t = v(:,1:r); % r-truncate U, Sigma & V
Kt = u_t*sig_t*v_t' ;   %Generates the truncated Kt
%January mod: take away GAK


%  %January mod: new K
% [U,S,V] = svd(Psi_mX / sqrt(M), 'econ');  % Psi_mX is (2D x M)
% singvals = diag(S);
% cfg.trunk = 1e-4;
% r = find(singvals >= cfg.trunk, 1, 'last');
% U = U(:,1:r);
% S = S(1:r,1:r);
% V = V(:,1:r);
% Kt2 = (S \ (U' * (Psi_mY/sqrt(M)) * V));    

%Up to here, we know kt dimentions are given by the dictionary size and kt2
%by the truncation, even when their rank is the same. We need to find
%(downstream) what operation (similarity transformation of just a lifting)




%Calculate and normalize eigenvectors: 
[XI,MU,W] = eig(Kt,"nobalance");       %Generate eigenvalues and left-right eigenvctors
[XIb,MUb,Wb] = eig(Kt); %balanced eigenvalues

inerp = sum(conj(W).*XI);                  %Calculate the w_k'xi_k products
%inerp(inerp< 0.00001) = 0;   %Remove neglible inner p. (0.00001 worked ok) 
Wn = W./inerp;    %define the scaled w_n left eigenvectors in the u. circle
scaledIP = sum(conj(Wn).*XI);              %Calculate the w_n'xi_k products
ip_angle = angle(scaledIP);     %phase of w_n'xi_k products in the c. plane
%an_cor = exp(1i * +ip_angle);                      %angle correction factor
an_cor = exp(1i * ip_angle);                      %angle correction factor
Wnn = Wn.*an_cor;                                 %normalized w eigevectors
%Wnn2 = V*Wnn; %January mod
ipfinal = sum(conj(Wnn).*XI); %this should be vector of ones if all went ok

%Calculate the Koopman modes v_i
V = (Wnn'*B).';    %eDMDSVD: We use conjugate transpose now...
%V = (Wnn2'*B).';    %January mode
V(isnan(V)) = 0; %remove all undefined entries.

%State reconstruction with discrete-time eigenvalues:
mu = repmat(diag(MU)', M, 1);
row_indices = (1:M)';  % Create a column vector of row indices
mu_p = mu .^ row_indices;  % Element-wise exponentiation
Psi0 = Psi_mX(1,:);
Phid = diag(Psi0*XI);  %State reconstruction: build the diagonal Phi matrix
Xrr = real((mu_p*Phid*V')'); %State recontruction: 
orN = size(X_0_rec,1); %Original number of channels
Xr = Xrr(1:orN,:); %Remove copies from the stacking 
Xr_col = Xr(:);
Xr_col_n = normalize(Xr_col, "range");
Xr_rec = reshape(Xr_col_n, orN, M);

%ERROR:
%rem=Xr_rec-X_0(:, 1:1597); %This needs to be generalized!!!
%E=norm(rem,"fro")/norm(X1,"fro");


%Graphs: Artificial dataset and eDMD reconstruction
figure;
subplot(2,1,1)
surf(X_0_rec, 'EdgeColor','none')
title('artificial dataset', siggma)
subplot(2,1,2)
surf(Xr_rec, 'EdgeColor','none')
title(sprintf('eDMD, period (\\omega = %.2f)', omega));   % Add omega value
title('eDMD reconstruction')


%Graphs: SVD singular values (and truncated sv)
figure;
subplot(1,2,1)
plot(sig)
title('SVD SIGMA');
subplot(1,2,2)
plot(sig_t)
title('truncated SVD sigma');

