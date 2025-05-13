clear all
%PARAMETERS:
dt = 0.004;                  % delta*t Considering a sampling rate of 250Hz
nstacks = 10 ;                                          %Stacking parameter
r = 11;                                              %SVD: rank truncation 
N = 70;                 %ARTIFICIAL DATASET: number of channels in the data
M = 1600;              %ARTIFICIAL DATASET: number of snapshots in the data
t = dt.*(1:M);                             %ARTIFICIAL DATASET: Time vector
omega = 17.5;                      % ARTIFICIAL DATASET: Frequency multiplier
tau = 1.2;                            % ARTIFICIAL DATASET: delay parameter
Np = 0.0000;                               % ARTIFICIAL DATASET: Noise parameter
MA = 3;                      %ARTIFICIAL DATASET: Moving average parameter

%ARTIFICIAL DATASET:
No = Np.*(randn(N,M));                                         %White noise
f0 = 1;          % start frequency
f1 = 0.01;         %end freq
%X_osc = cos(2*pi*(f0*t + (f1 - f0)/(2*t(end)) * t.^2));
%X_osc = t+cos(2*pi*t);
%X_osc = t.^2;
%X_osc = 1./(1+exp(-t));
X_osc = t;
X_zero= zeros(N,M);
X_zero(5:7, :) = repmat(0.7.*X_osc, 3, 1);
X_zero(17:19, :) = repmat(0.95.*X_osc, 3, 1);
X_zero(55:57, :) = repmat(-0.95.*X_osc, 3, 1);
X_zero(65:67, :) = repmat(-0.35.*X_osc, 3, 1);
%Xconst= zeros(N,M);
%Xconst(5:7, :) = 2;
X_0 = X_zero + No;
Xcol = X_0(:);
X_col_n = normalize(Xcol, "range");
X_0_rec = reshape(X_col_n, N, M);
%X_0 = movmean(X_0,MA);   %This is a mobile mean function that makes our data a bit more smoth
X1 = X_0_rec(:, 1:end-1);
X2 = X_0_rec(:, 2:end);


%Here we will perform stacking:
%nstacks = 1;  %Stacking, at the parameters section
if nstacks > 1
    Xaug = [];
    for st = 1:nstacks
        Xaug = [Xaug; X_0(:, st:end-nstacks+st)]; 
    end
    
    X1 = Xaug(:, 1:end-1);
    X2 = Xaug(:, 2:end);
else
    X1 = X_0(:, 1:end-1);
    X2 = X_0(:, 2:end);
end
N = size(X1,1);  %number of channels in the stack
M = size(X1,2); %number of snapshots in the stack 

%A very very short eDMD dictionary, for testing purposes
% X1e = [0.*abs(X1).^2; X1];
% X2e = [0.*abs(X2).^2; X2];

%We define activity clusters to help us represent activity relations
C1 = mean(X_0(5:7, 1:end-nstacks), 1);
C2 = mean(X_0(17:19, 1:end-nstacks), 1);
C3 = mean(X_0(55:57, 1:end-nstacks), 1); %unrelated area in the occipital lobe


% X1e = [X1 ; 0.*C1.*X1 ; 0.*C2.*X1 ; 0.*C3.*X1 ];
% X2e = [X2 ; 0.*C1.*X2 ; 0.*C1.*X2 ; 0.*C1.*X2 ];
% X1e = [X1 ; 0.*C1.*X1 ; 0.*C2.*X1 ; 0.*C3.*X1 ; 0.*asind(X1); 0.*acos(X1) ];
% X2e = [X2 ; 0.*C1.*X2 ; 0.*C1.*X2 ; 0.*C1.*X2 ; 0.*asind(X2); 0.*acos(X2) ];
% X1e = [X1 ; C1.*X1 ; C2.*X1 ; C3.*X1 ; asind(X1); acos(X1) ];
% X2e = [X2 ; C1.*X2 ; C1.*X2 ; C1.*X2 ; asind(X2); acos(X2) ];
%X1e = [X1 ; 0.* cos(2*pi*(f0*X1 + (f1 - f0)/(2*X1(end)) * X1.^2)) ];
%X2e = [X2 ; 0.* cos(2*pi*(f0*X2 + (f1 - f0)/(2*X2(end)) * X2.^2))  ];
% X1e = [X1 ;  cos(2*pi*(f0*X1 + (f1 - f0) * X1.^2)) ];
% X2e = [X2 ;  cos(2*pi*(f0*X2 + (f1 - f0) * X2.^2)) ];
% X1e = [X1 ; 1./(1+exp(-X1))  ];
% X2e = [X2 ; 1./(1+exp(-X2))  ];

X1e = [X1 ];
X2e = [X2 ];



Psi_mX = X1e.'; %Proper variable naming
Psi_mY = X2e.';

%THE B MATRIX: weights that are used to build the koopman modes
%It is easily built by realizing X1 in X1e encondes the full state
k = size(X2e, 1);  %dimension of the observables vector
% btop = zeros(N,N);
% bmid = eye(N);
% B = [btop;bmid];


%b_top = eye(N);
%b_low = zeros(1*N,N);
%B = [b_top;b_low];
B = eye(N);
%X1e and X2e are Psi_M transposed, run the next couple lines to verify
%that B fulfils eq16 in Williams:
 %  Cr=transpose(B)*transpose(Psi_mX);
 %  CCr=Cr-X1;


%X1e and X2e are Psi_M transposed, run the next couple lines to verify
%that B fulfils eq16 in Williams:
%  C=transpose(B)*transpose(Psi_mX);
%  CC=C-X1;

%Generate the eDMD approximation to the Koopman Operator:
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
[XI,MU,W] = eig(Kt);       %Generate eigenvalues and left-right eigenvctors
inerp = sum(conj(W).*XI);                  %Calculate the w_k'xi_k products
%inerp(inerp< 0.00001) = 0;   %Remove neglible inner p. (0.00001 worked ok) 
Wn = W./inerp;    %define the scaled w_n left eigenvectors in the u. circle
scaledIP = sum(conj(Wn).*XI);              %Calculate the w_n'xi_k products
ip_angle = angle(scaledIP);     %phase of w_n'xi_k products in the c. plane
an_cor = exp(1i * +ip_angle);                      %angle correction factor
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
orN = size(X_0,1); %Original number of channels
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
