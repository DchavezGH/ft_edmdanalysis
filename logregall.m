theta_pvalue_session = [];    %stored p_values 
subjects= 1:24;
subjects(subjects == 13) = [];

for i = subjects


freqmat_labeled_s = [];
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
freqmat_labeled_s = [freqmat_labeled_s; freqmat_labeled];

end  %if
end  %if
end  %for: j 1 to 4


T = array2table(freqmat_labeled,'VariableNames',{'Delta','Theta','Alpha','lowBeta','highBeta','Gamma','responsevar'});
modelspec = "responsevar ~   Theta*lowBeta*highBeta";  
mdl = fitglm(T,modelspec, Distribution="binomial") %LOG REG happens here
LRcoeffs = mdl.Coefficients;
%session_theta_pvalue(i,:) =LRcoeffs{2,4} ;    %stored p_values 
%session_l_beta_pvalue(i,:) =LRcoeffs{2,4} ;    %stored p_values 
%session_h_beta_pvalue(i,:) =LRcoeffs{2,4} ;    %stored p_values 
session_mixed_pvalue(i,:) =LRcoeffs{6,4} ;    %stored p_values


end 

% PLOTS % PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS
% PLOTS % PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS
% PLOTS % PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS% PLOTS

p_vals=[session_theta_pvalue session_l_beta_pvalue session_h_beta_pvalue session_mixed_pvalue];
p_vals(13,:)=[];
p_vals=p_vals';


figure;

% Plot with ball markers
%h = plot(p_vals, 'LineWidth', 0.75, ...
%         'Marker', 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'auto');
h = plot(p_vals, ...
    'LineStyle', 'none', ...       % ← removes the connecting lines
    'Marker', 'o', ...
    'MarkerSize', 4, ...
    'MarkerFaceColor', 'auto');



%legend({'Theta', 'Low Beta', 'High Beta'}, ...
%       'Location', 'northoutside', 'Box', 'off', ...
%      'Orientation', 'horizontal', IconColumnWidth=10);

%xlabel('Band', 'FontName', 'Helvetica', 'FontSize', 10);
ylabel('pValue per subject', 'FontName', 'Helvetica', 'FontSize', 10);

% Axis formatting
ax = gca;
ax.FontName = 'Helvetica';
ax.FontSize = 10;
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'Theta', 'Low Beta', 'High Beta', 'Theta-Beta'};
xlim([0.5 4.5]);  % now columns 1–4 are visible

%xlim([0.5 3.5]);   % <-- adds spacing on both sides
ax.YTick = [ 0 .25    1];
ylim padded 
box on

% Figure size
width = 76/10;   % mm → cm
height = 50/10;
set(gcf, 'Units', 'centimeters', 'Position', [5 5 width height]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [width height]);

% Grey dotted baseline
yline(0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
yline(0.25, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
% Export as vector PDF
print(gcf, 'pvalsmx.pdf', '-dpdf', '-painters');



