clear;
close all;
thisDir = fileparts(mfilename('fullpath'));

budgetList = [5,10,15,20,25];      % default: 10
windowList = [5,10,15,20,25];      % default: 15
sparseList = [1,2,3,4,5];          % default: 1
noiseList  = [0.01,0.05,0.1,0.5,1]; % default: 0.01

meanBudget = zeros(numel(budgetList),1);
stdBudget  = zeros(numel(budgetList),1);
meanWindow = zeros(numel(windowList),1);
stdWindow  = zeros(numel(windowList),1);
meanSparse = zeros(numel(sparseList),1);
stdSparse  = zeros(numel(sparseList),1);
meanNoise  = zeros(numel(noiseList),1);
stdNoise   = zeros(numel(noiseList),1);

for j = 1:numel(budgetList)
    S = load([thisDir,'\0ARL_save\Sensitivity_budget_m\SB-FQCD_','budget',num2str(j)]);
    meanBudget(j) = S.meanDelay;
    stdBudget(j)  = S.stdDelay;
end

for j = 1:numel(windowList)
    S = load([thisDir,'\0ARL_save\Sensitivity_window_w\SB-FQCD_','window',num2str(j)]);
    meanWindow(j) = S.meanDelay;
    stdWindow(j)  = S.stdDelay;
end

for j = 1:numel(sparseList)
    S = load([thisDir,'\0ARL_save\Sensitivity_sparsity_s\SB-FQCD_','sparse',num2str(j)]);
    meanSparse(j) = S.meanDelay;
    stdSparse(j)  = S.stdDelay;
end

for j = 1:numel(noiseList)
    S = load([thisDir,'\0ARL_save\Sensitivity_noise_sigma\SB-FQCD_','noise',num2str(j)]);
    meanNoise(j) = S.meanDelay;
    stdNoise(j)  = S.stdDelay;
end

%% Unified style (Scheme A)
figPos = [100 100 260 210];
baseDir = fullfile(fileparts(thisDir), 'figures');
if ~exist(baseDir, 'dir')
    mkdir(baseDir);
end

lineColor = [0 0 0];           % black line
shadeColor = [0.75 0.75 0.75]; % gray shaded area
markerType = 'o';
lineWidth = 1.5;
markerSize = 3;
shadeAlpha = 0.35;

%% 1. Budget (m)
figure;
set(gcf, 'Position', figPos);
x = 1:numel(budgetList);
plot_with_shade(x, meanBudget, stdBudget, shadeColor, shadeAlpha, ...
    lineColor, markerType, lineWidth, markerSize);
xlim([x(1)-0.2, x(end)+0.2]);
xticks(x);
xticklabels(string(budgetList));
ylabel('ADD');

box on;


fullPath = fullfile(baseDir, 'add_budget.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);

%% 2. Window (w)
figure;
set(gcf, 'Position', figPos);
x = 1:numel(windowList);
plot_with_shade(x, meanWindow, stdWindow, shadeColor, shadeAlpha, ...
    lineColor, markerType, lineWidth, markerSize);
xlim([x(1)-0.2, x(end)+0.2]);
xticks(x);
xticklabels(string(windowList));
ylabel('ADD');
box on;

fullPath = fullfile(baseDir, 'add_window.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);

%% 3. Sparse (s)
figure;
set(gcf, 'Position', figPos);
x = 1:numel(sparseList);
plot_with_shade(x, meanSparse, stdSparse, shadeColor, shadeAlpha, ...
    lineColor, markerType, lineWidth, markerSize);
xlim([x(1)-0.2, x(end)+0.2]);
xticks(x);
xticklabels(string(sparseList));
ylabel('ADD');
box on;


fullPath = fullfile(baseDir, 'add_sparse.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);

%% 4. Noise (\sigma^2)
figure;
set(gcf, 'Position', figPos);
x = 1:numel(noiseList);
plot_with_shade(x, meanNoise, stdNoise, shadeColor, shadeAlpha, ...
    lineColor, markerType, lineWidth, markerSize);
xlim([x(1)-0.2, x(end)+0.2]);
xticks(x);
xticklabels(string(noiseList));
ylabel('ADD');
box on;


fullPath = fullfile(baseDir, 'add_noise.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);

%% Local function
function plot_with_shade(x, meanVal, stdVal, shadeColor, shadeAlpha, ...
    lineColor, markerType, lineWidth, markerSize)

    hold on;

    upper = meanVal(:)' + stdVal(:)';
    lower = meanVal(:)' - stdVal(:)';

    % shaded band: mean +/- std
    fill([x, fliplr(x)], [upper, fliplr(lower)], shadeColor, ...
        'FaceAlpha', shadeAlpha, 'EdgeColor', 'none');

    % main line
    plot(x, meanVal, ['-' markerType], ...
        'Color', lineColor, ...
        'MarkerFaceColor', lineColor, ...
        'MarkerEdgeColor', lineColor, ...
        'LineWidth', lineWidth, ...
        'MarkerSize', markerSize);

    hold off;
end





