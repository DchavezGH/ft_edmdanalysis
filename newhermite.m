clear all
%PARAMETERS:
%dt = 0.004;                  % delta*t Considering a sampling rate of 250Hz
nstacks = 5;                                          %Stacking parameter
r = 30;                                              %SVD: rank truncation 
N = 8;                 %ARTIFICIAL DATASET: number of channels in the data
M = 160;              %ARTIFICIAL DATASET: number of snapshots in the data
%t = dt.*(1:M);       %ARTIFICIAL DATASET: Time vector, temporarily ignored
t = linspace(0, 1, M); %This is needed further down to build B
omega = 17.5;                      % ARTIFICIAL DATASET: Frequency multiplier
tau = 1.2;                            % ARTIFICIAL DATASET: delay parameter
MA = 3;                      %ARTIFICIAL DATASET: Moving average parameter

%ARTIFICIAL DATASET:
Nn= N;     %Safekeeping original channel dim
Xs = 0.9.*exp(-(20.*t-1.8).^2);
Xss = -0.9.*exp(-(20.*t-0.8).^2);
Xi = repmat(Xs, N/2, 1);
Xii = repmat(Xss, N/2, 1);
X_zero = [Xi;Xii];
X_0_rec = X_zero; %take it off when retaking normalization




%We temporalily remove normalization
% X_0 = X_zero + No;
% Xcol = X_0(:);
% X_col_n = normalize(Xcol, "range");
% X_0_rec = reshape(X_col_n, N, M);
% %X_0 = movmean(X_0,MA);   %This is a mobile mean function that makes our data a bit more smoth
% X1 = X_0_rec(:, 1:end-1);
% X2 = X_0_rec(:, 2:end);




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

X1e = [ X1.^2-1 ; X1.^3-3.*X1 ; X1.^4-6.*X1.^2+3 ; X1.^5-10.*X1.^3+15.*X1 ; X1.^6-15.*X1.^4+45.*X1.^2-15 ; X1.^7-21.*X1.^5+1.5.*X1.^3-105.*X1 ; X1.^8-28.*X1.^6+210.*X1.^4-420.*X1.^2+105 ; X1.^9-36.*X1.^7+378.*X1.^5-1260.*X1.^3+495.*X1];
X2e = [ X2.^2-1 ; X2.^3-3.*X2 ; X2.^4-6.*X2.^2+3 ; X2.^5-10.*X2.^3+15.*X2 ; X2.^6-15.*X2.^4+45.*X2.^2-15 ; X2.^7-21.*X2.^5+1.5.*X2.^3-105.*X2 ; X2.^8-28.*X2.^6+210.*X2.^4-420.*X2.^2+105 ; X2.^9-36.*X2.^7+378.*X2.^5-1260.*X2.^3+495.*X2]; 

Psi_mX = X1e.'; %Proper variable naming
Psi_mY = X2e.';

%GENERATE THE B MATRIX: weights that are used to build the koopman modes

degree = 9; %Dictionary degree
%t_stack = dt.*(1:M);   %t vector for the dims of the stack

% Build design matrix (T x (degree+1)) for the time domain
Phi = Psi_mX;

%Phi should actually contain the elements of Psi_mY, so the least squares
%problem delivers a useful coeff. matrix: 



% Preallocate coefficients
A = zeros(N, (degree-1)*(N));   % Rows: channels, Columns: polynomial terms

% Fit each channel's time-series data using least squares
for m = 1:N
    y = X1(m, :)';                 % Column vector for time series of channel m
    a = Phi \ y;                  % Solve least squares: Phi * a ≈ y
    A(m, :) = a';                 % Store in row
end




lqrec = Phi * A';

B=pinv(Phi)*X1';
g=B.'*Psi_mX.';


%Run the next couple lines to verify
%that B fulfils eq16 in Williams:
g=B.'*Psi_mX.';   %the full state observable g(x)
%CC=g-X1;


figure;
surf(g, 'EdgeColor','none')
title('Full state observabe g(x)')






G = (Psi_mX'*Psi_mX)/M;
A = (Psi_mX'*Psi_mY)/M;
K = pinv(G)*A;    %The eDMD approximation to the Koopman operator
clear G A

%We introduce SVD to truncate K
[u,sig,v] = svd(K, 'econ'); % Performs svd on the data
%r=5; %rank truncation parameter is at the parameters section
u_t = u(:,1:r); 
sig_t = sig(1:r,1:r); 
v_t = v(:,1:r); % r-truncate U, Sigma & V
Kt = u_t*sig_t*v_t' ;   %Generates the truncated Kt

%Calculate the and normalize eigenvectors: 
[XI,MU,W] = eig(Kt,"nobalance");       %Generate eigenvalues and left-right eigenvctors
[XIb,MUb,Wb] = eig(Kt); %balanced eigenvalues

%balanced/unbalanced discrepancies
XIr_d = real(XI)-real(XIb); 
XIi_d = imag(XI)-imag(XIb);
MUr_d = real(MU)-real(MUb); 
MUi_d = imag(MU)-imag(MUb);
Wr_d = real(W)-real(Wb);
Wi_d = imag(W)-imag(Wb);

E_XIr=norm(XIr_d ,"fro")/norm(real(XI),"fro");
E_XIi=norm(XIi_d ,"fro")/norm(imag(XI),"fro");

E_MUr=norm(MUr_d ,"fro")/norm(real(MU),"fro");
E_MUi=norm(MUi_d ,"fro")/norm(imag(MU),"fro");

E_Wr=norm(Wr_d ,"fro")/norm(real(W),"fro");
E_Wi=norm(Wi_d  ,"fro")/norm(imag(W),"fro");

inerp = sum(conj(W).*XI);                  %Calculate the w_k'xi_k products
%inerp(inerp< 0.00001) = 0;   %Remove neglible inner p. (0.00001 worked ok) 
Wn = W./inerp;    %define the scaled w_n left eigenvectors in the u. circle
scaledIP = sum(conj(Wn).*XI);              %Calculate the w_n'xi_k products
ip_angle = angle(scaledIP);     %phase of w_n'xi_k products in the c. plane
%an_cor = exp(1i * +ip_angle);                      %angle correction factor
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
title('artificial dataset')
subplot(2,1,2)
surf(Xr_rec, 'EdgeColor','none')
title(sprintf('eDMD, period (\\omega = %.2f)', omega));   % Add omega value
title('eDMD reconstruction')

%Graphs: Eigenvalues, left and right eigenvectors
figure;
subplot(2,3,1)
surf(real(Wnn), 'EdgeColor','none')
title('real W normalized');
subplot(2,3,2)
surf(real(XI), 'EdgeColor','none')
title('real XI');
subplot(2,3,3)
surf(real(MU), 'EdgeColor','none')
title('real MU');
subplot(2,3,4)
surf(imag(Wnn), 'EdgeColor','none')
title('imag W normalized');
subplot(2,3,5)
surf(imag(XI), 'EdgeColor','none')
title('imag XI');
subplot(2,3,6)
surf(imag(MU), 'EdgeColor','none')
title('imag MU');

%Graphs: SVD singular values (and truncated sv)
figure;
subplot(1,2,1)
plot(sig)
title('SVD SIGMA');
subplot(1,2,2)
plot(sig_t)
title('truncated SVD sigma');
