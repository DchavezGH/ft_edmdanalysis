%This code performs extendedDMD on a EEG dataset and calculates the error
%between the generated eDMD representation and the original data.
%We follow the eDMD method of Williams, Kevrekidis and Rowley (2014), but
%the algorithm implementation is our own. 
%To read the dataset, we used a bit of code from Makoto (https://sccn.ucsd.edu/wiki/Makoto%27s_useful_EEGLAB_code)
%The dataset is taken from: 
%   (https://openneuro.org/datasets/ds004117/versions/1.0.1/file-display/participants.json)
%and reflects a WM paradigm by Julie Onton's team 



%LOAD DATASET:
%Load header info:
headerInfo = load('sub-001_ses-01_task-WorkingMemory_run-1_eeg.set', '-mat');
% Load EEG time-series matrix (taken from eeg_getdatact.m centered at line 233)
fid = fopen('sub-001_ses-01_task-WorkingMemory_run-1_eeg.fdt', 'r', 'ieee-le');
for trialIdx = 1:headerInfo.trials % In case the saved data are epoched, loop the process for each epoch. Thanks Ramesh Srinivasan!
    currentTrialData = fread(fid, [headerInfo.nbchan headerInfo.pnts], 'float32');
    data(:,:,trialIdx) = currentTrialData; % Data dimentions are: electrodes, time points, and trials (the last one is for epoched data)                  
end
fclose(fid);
X_0 = currentTrialData; %Raw data
X_0 = movmean(X_0,10);     %This is a mobile mean function that makes our data a bit more smoth

%SET DMD PARAMETERS:
dt = 0.004; % delta*t Considering a sampling rate of 250Hz
%dt = 0.4; % slower sampling for testing purposes
N = size(X_0, 1);  %number of channels in the data
M = fix(size(X_0, 2)); %number of snapshots in the data
t = dt.*(1:M); %Time vector

%ARTIFICIAL DATASET:
X_osc = 10.* sin(repmat(0.1.*t, 71, 1)); %This is a sine oscillation, comment it when not needed
Xconst = zeros(N,M);
    %Xconst (1:10, :) = 10;
    %X_0=X_0+X_osc;   %Here, we instroduce the sine oscillation that should be trivial to represent by DMD
X_0=X_osc+Xconst;

X1 = X_0(:, 1:end-1);
X2 = X_0(:, 2:end);

%A very very short eDMD dictionary, for testing purposes
%X1e=[X1; abs(X1).^2; ];
%X2e=[X2; abs(X2).^2];


%A eDMD dictionary made of the 10 first Hermite Polynomials.
%MATLAB's Hermite function is symbolic, we recommend against using it.
%Strongly consider building a dictionary with the products of 
%Hermite Polynomials, but the computational cost will grow a lot
%We stop using it since the sheer magnitude of the last hermite polynomials
%make it very hard to operate.
H0 = ones(size(X1));                           
H1 = 2 * X1;   
H2 = 4 * X1.^2 - 2;                            
H3 = 8 * X1.^3 - 12 * X1;                      
H4 = 16 * X1.^4 - 48 * X1.^2 + 12;             
H5 = 32 * X1.^5 - 160 * X1.^3 + 120 * X1;      
H6 = 64 * X1.^6 - 480 * X1.^4 + 720 * X1.^2 - 120; 
H7 = 128 * X1.^7 - 1344 * X1.^5 + 3360 * X1.^3 - 1680 * X1; 
H8 = 256 * X1.^8 - 3584 * X1.^6 + 13440 * X1.^4 - 13440 * X1.^2 + 1680; 
H9 = 512 * X1.^9 - 9216 * X1.^7 + 48384 * X1.^5 - 80640 * X1.^3 + 30240 * X1; 
H10 = 1024 * X1.^10 - 23040 * X1.^8 + 161280 * X1.^6 - 403200 * X1.^4 + 302400 * X1.^2 - 30240; 
X1e=[H0; H1; H2; H3; H4; H5; H6; H7; H8; H9; H10];

H0 = ones(size(X2));                           
H1 = 2 * X2; 
H2 = 4 * X2.^2 - 2;                            
H3 = 8 * X2.^3 - 12 * X2;                      
H4 = 16 * X2.^4 - 48 * X2.^2 + 12;             
H5 = 32 * X2.^5 - 160 * X2.^3 + 120 * X2;      
H6 = 64 * X2.^6 - 480 * X2.^4 + 720 * X2.^2 - 120; 
H7 = 128 * X2.^7 - 1344 * X2.^5 + 3360 * X2.^3 - 1680 * X2; 
H8 = 256 * X2.^8 - 3584 * X2.^6 + 13440 * X2.^4 - 13440 * X2.^2 + 1680; 
H9 = 512 * X2.^9 - 9216 * X2.^7 + 48384 * X2.^5 - 80640 * X2.^3 + 30240 * X2; 
H10 = 1024 * X2.^10 - 23040 * X2.^8 + 161280 * X2.^6 - 403200 * X2.^4 + 302400 * X2.^2 - 30240; 
X2e=[H0; H1; H2; H3; H4; H5; H6; H7; H8; H9; H10 ];

%THE B MATRIX: weights that are used to build the koopman modes
%It is easily built by realizing H1 in X2e encondes the full state
k = size(X2e, 1);  %dimention of the observables vector
btop = zeros(N,N);
bmid= eye(N)./2;
blow = zeros(k-2*N,N);
B=[btop;bmid;blow];
%X1e and X2e are Psi_M transposed, run the next couple lines to verify
%that B fulfils eq14 in Williams:
% C=transpose(B)*X1e;
% CC=C-X1;

%eDMD approximation to the Koopman Operator "AY" as in Kutz:
%A1= X2e*X1e.';
%A2= X1e*X1e.';
%AY=A1*pinv(A2);
%This algorithm did not work. We decide to implement instead William's 

%Generate the eDMD approximation to the Koopman Operator K as in Williams:
%eDMDSVD: We decide to transpose first X1e and X2e to avoid confusions:
%eDMDSVD: We use conjugate transpose ' instead of .'
Psi_mX = X1e';
Psi_mY = X2e';
G=(Psi_mX'*Psi_mX)/M;
A=(Psi_mX'*Psi_mY)/M;
K=pinv(G)*A;    %The eDMD approximation to the Koopman operator
% %%eDMDSVD: We introduce SVD to truncate K and hopefully fix the algorithm...
%   [u,sig,v] = svd(K, 'econ'); %%eDMDSVD: Performs svd on the data
%   r=40; %% rank truncation parameter
%   u_t=u(:,1:r); sig_t=sig(1:r,1:r); v_t=v(:,1:r); %eDMDSVD:  r-truncate U, Sigma & V
%   Kt =u_t*sig_t*v_t' ;   %eDMDSVD: generates the truncated Kt, but it's not
% %that helpful
% %Aqui el viernes pensamos en introducir la K truncada: 
% %Instead of Kt, define a new K=U*AU=U*Psi_mY*V*sig^{-1}
% Kt = u_t'*Psi_mY'*v_t/sig_t;


[XI,MU,W] = eig(K);   %eDMDSVD: Calculate (L and R) eigenvectors and eigenvalues 
W=inv(XI)';        %as an alternative to fix V, let's calculate W as in W19
  
XI(XI < 0.000001) = 0; %Let every neglible element in XI be zero
W(W < 0.000001) = 0; %Let every neglible element in XI be zero

%Normalize the left eigenvectors:
inerp=sum(conj(W).*XI);  %Calculate inner product across columns
%rinerp=real(inerp); %given how neglible is the imaginary part, we remove it
                     %we don't do that anymore since the imaginary part is
                     %not neccesarily neglible
inerp(inerp < 0.000001) = 0; %Let every neglible element in rinerp be zero
Wn = W./inerp;        % then divide w_i by such inner products to normalize them
inerpchk=sum(conj(Wn).*XI); %inerpchk should be a vector of 1s if all went ok. 


 figure;
plot(real(inerpchk));
 figure;
plot(imag(inerpchk));
figure;
plot(inerpchk);
                      %Working now. Naturally, entries involving zero
                      %eigenvalues are undefined
%Calculate the Koopman modes v_i
V=(Wn'*B).';    %eDMDSVD: We use conjugate transpose now...
        %%%%%%%AUDITED UP TO HERE
V(isnan(V)) = 0;
 
 figure;
  surf(Psi_mX, 'EdgeColor','none')  
  title('W')        
figure;
  surf(imag(XI), 'EdgeColor','none')  
  title('XI')



%State reconstruction: time vector t is now a diagonal matrix: 
% do not use diag(t), as it takes too much memory
diat = spdiags(t(:), 0, length(t), length(t));
%State reconstruction: build the elt:[e^lambda*t] matrix: 
mu = repmat(diag(MU)', M, 1);
elt = exp(diat*(log(mu)./dt));                              %exp^Lambda*t
%State reconstruction: build the diagonal Phi matrix
Psi0 = X1e(:,1)';
Phid = diag(Psi0*XI);                                       %Phi diagonal
                                               %V matrix
%State recontruction: 
Xr = real((elt*Phid*V')');
  figure;
  surf(real(Xr), 'EdgeColor','none')  
  title('reconstruction')



%ERROR(has to be fixed with the new variables):
%rem=X1-Xr;
%E=norm(rem,"fro")/norm(X1,"fro")

%GRAPHS (has to be fixed with the new variables):
figure;
imagesc(real(Xr));
title('Xr')
figure;
imagesc(X1);
title('rawData')
% figure;
% imagesc(rem);
% title('error')

%3D GRAPHS, THEY TEND TO BE CUMBERSOME BUT SOMETIMES INFORMATIVE(have to be fixed with the new variables):
%  figure;
%  surf(X_0, 'EdgeColor','none')
%  title('eDMD')
% figure;
% surf(X1, 'EdgeColor','none')
% title('rawData')
%figure;
%surf(rem, 'EdgeColor','none')
% title('ERROR')

%figure;
%imagesc(Kt(:,1));

  figure;
  surf(imag(V), 'EdgeColor','none')  
  title('Vr')
  figure;
  surf(K, 'EdgeColor','none')