%This code:
%  -Requires the total power per band generated from alldata_eDMD.m
%  -Takes the theta power per session (asumes a given sessions keeps the
%  same EEG configuration), locates the 5-letter trials corresponding to
%  the [90 70 50 30]th percentiles and generates Pearson Correlation
%  coefficients between theta_power and number of encoded letters

for kk = [1 2 3 4]  %used percentiles
tille= [90 70 50 30];
subjects= 1:24;
subjects(subjects == 13) = [];

for i = subjects
theta_cor_session = [];    
theta_power_session =[];    

for j = 1:4
if isempty(peakfeats{i, j}) %Get rid of empty cells
sumfeatures = NaN; 
else
                   
sumfeatures = peakfeats{i, j};

vl = size(sumfeatures, 1); %how many events
n = vl/5;                  %how many trials
theta_cor =zeros(vl,2);

sumfeatures(:, [1]) = [];        %Remove all but theta power values
sumfeatures(:, [2 3 4 5]) = NaN;  %Remove all but theta power values
sumfeatures(:,2) = repmat((1:5).', n, 1); %2nd col: how many encoded letters
%Calculate per trial power
for ii = 1:5:size(sumfeatures,1)
    idx = ii:ii+4;                 % indices for this block of 5
    s = sum(sumfeatures(idx,1));           % sum of 5 elements in column 1
    sumfeatures(idx,3) = s;                % store the sum in column 3
end

 %agregate the whole sesion in an array

theta_cor_session = [theta_cor_session; sumfeatures];


end %if 
end %for j 1 to 4

%Calculate kk-th percentile
%threshold = prctile(theta_cor_session(:,3), 90); % Compute the 90th percentile threshold
threshold = prctile(theta_cor_session(:,3), tille(kk)); % Compute the 90th percentile threshold
theta_cor_session(:,4) = theta_cor_session(:,3) >= threshold;
theta_cor_session = theta_cor_session(theta_cor_session(:,4) == 1, :); %90th percentile trials of the session, together
%theta_cor_session(:,4) = theta_cor_session(theta_cor_session >= threshold); % Extract elements in the 90th percentile
R = corrcoef(theta_cor_session);
% local_corr(i) = R(1,2);
local_corr(i,kk) = R(1,2);

end 
end
 
 local_corr(13,:) = []; %as this one is empty

%clearvars -except peakfeats



%PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT 
figure;

% Plot with ball markers
h = plot(local_corr, 'LineWidth', 0.75, ...
         'Marker', 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'auto');

legend({'90th', '70th', '50th', '30th'}, ...
       'Location', 'northoutside', 'Box', 'off', ...
       'Orientation', 'horizontal', IconColumnWidth=10);

xlabel('Subject', 'FontName', 'Helvetica', 'FontSize', 10);
ylabel('Correlation Coefficient', 'FontName', 'Helvetica', 'FontSize', 10);

% Axis formatting
ax = gca;
ax.FontName = 'Helvetica';
ax.FontSize = 10;
ax.XTick = [1 5 10 15 20];
ax.YTick = [-.5 -.25 0 .25 .5];
box on

% Figure size
width = 76/10;   % mm → cm
height = 50/10;
set(gcf, 'Units', 'centimeters', 'Position', [5 5 width height]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [width height]);

% Grey dotted baseline
yline(0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
yline(0.25, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
yline(-0.25, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
% Export as vector PDF
print(gcf, 'pearsons2.pdf', '-dpdf', '-painters');