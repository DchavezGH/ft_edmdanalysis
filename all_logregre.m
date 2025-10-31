% %Artificial features
% tbl = zeros(10,9);
% tbl(:,1:8) = features;
% tbl(6:10,9)= 1;
% T = array2table(tbl,'VariableNames',{'Delta','Theta','Alpha','lowBeta','highBeta','Gamma','responsevar'});

%sumfeats_control=sumfeats; %run this once to get the control events apart


%for i = subjects
for i = 1            %remove
freqmat_labeled_s = [];
theta_pvalue_session = [];    %stored p_values 
highbeta_pvalue_session = [];    %stored p_values 

for j = 1:4
freqmat = [];    %sesion/subject dominant freq array
if isempty(sumfeats{i, j}) %Get rid of empty cells
freqmat = NaN;
else

if isempty(sumfeats_control{i, j}) %Get rid of empty cells
freqmat_c = NaN;
else          

lr_features = sumfeats{i, j}; %grab the j-th encoded run cell
lr_features_c = sumfeats_control{i, j}; %grab the j-th control run cell

% lr_features = sumfeats{1, 1}; %remove
% lr_features_c = sumfeats_control{1, 1}; %remove


freqmat = [lr_features; lr_features_c];
labels = [ones(length(lr_features), 1); zeros(length(lr_features_c), 1)];
freqmat_labeled= [freqmat, labels]; %1 means encoding event, 0 means control event
freqmat_labeled_s = [freqmat_labeled_s,freqmat_labeled];

end  %if
end  %if
end  %for: j 1 to 4

T = array2table(freqmat_labeled,'VariableNames',{'Delta','Theta','Alpha','lowBeta','highBeta','Gamma','responsevar'});
modelspec = "responsevar ~   highBeta";  
mdl = fitglm(T,modelspec, Distribution="binomial") %LOG REG happens here
theta_pvalue_session(i,j) =LRcoeffs{2,4} ;    %stored p_values 

end  %for: subjects 1 to 24







tbl = zeros(30,7);
tbl(:,1:6) = peakfeatures(:,1:6);
tbl(1:15,7)= 1;  

%T = array2table(tbl,'VariableNames',{'alpha','responsevar'});
T = array2table(tbl,'VariableNames',{'Delta','Theta','Alpha','lowBeta','highBeta','Gamma','responsevar'});
modelspec = "responsevar ~   Theta";  
mdl = fitglm(T,modelspec, Distribution="binomial");
LRcoeffs = mdl.Coefficients;
theta_pvalue_session(i,j) =LRcoeffs{2,4} ;    %stored p_values 

% 
% 
% 
% %PLOTS are still to be modified: 
% idxSmoker = T.responsevar == 1;
% idxNonSmoker = T.responsevar == 0;
% 
% 
% figure;
% hold on;
% 
% % Non-smokers: green crosses, larger and thicker
% plot(T.alpha(idxNonSmoker), T.beta(idxNonSmoker), ...
%     'gx', 'MarkerSize', 10, 'LineWidth', 1.5);
% 
% % Smokers: red circles, larger and thicker
% plot(T.alpha(idxSmoker), T.beta(idxSmoker), ...
%     'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
% 
% xlabel('alpha');
% ylabel('beta');
% title('Systolic vs Diastolic by Smoker Status');
% legend({'Non-smoker','Smoker'});
% grid on;
% hold off;
% 
% plotSlice(mdl)
% 



% figure;
% hold on;
% 
% % Non-smokers: green crosses, larger and thicker
% plot(T.alpha(idxNonSmoker), T.beta(idxNonSmoker), ...
%     'gx', 'MarkerSize', 10, 'LineWidth', 1.5);
% 
% % Smokers: red circles, larger and thicker
% plot(T.alpha(idxSmoker), T.beta(idxSmoker), ...
%     'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
% 
% xlabel('alpha');
% ylabel('beta');
% title('Systolic vs Diastolic by Smoker Status');
% legend({'Non-smoker','Smoker'});
% grid on;
% hold off;

