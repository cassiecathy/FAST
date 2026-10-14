clc;
clear;
close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(thisDir));

sigma2 = 0.01;
tau0 = 50;
T = 100;
w = 15;
initialRandomSteps = 20;
shift = 5;
thresholdPlot = inf;

M = 100;
d = 20;
m = 10;
Sigma = eye(d);
Phi = build_basis_matrix(M, d);
betaSupport = [10];
betaSigns = [1];

rng('default');

%% time
runtimeMat = nan(T, 3);
WPath_FUCB = nan(T, 3);
betaErr_FUCB = nan(T, 3);
regret_FUCB = nan(T, 3);

scaleM = [20, 100, 500];
scaleD = [8, 20, 50];
scalem = [4,10,25];
scalew = [10,15,20];
betaSupportk = [3,10,25];
betaSignsk = 1;

for k = 1:3
    Mk = scaleM(k);
    dk = scaleD(k);
    mk = scalem(k);
    wk = scalew(k);
    Sigmak = eye(dk);
    Phik = build_basis_matrix(Mk, dk);

    [~, WPath_FUCB(:,k), betaErr_FUCB(:,k), regret_FUCB(:,k), ~, runtimeMat(:, k)] = sparse_bandit( ...
        Phik, Sigmak, sigma2, shift, tau0, T, mk, wk, inf, ...
        'proposed', initialRandomSteps, 'beta', betaSupportk(k), betaSignsk);
end

%% plot figures separately
betaErr_FUCB(1:tau0,:) = [];
regret_FUCB(1:tau0,:) = [];
WPath_FUCB(1:tau0,:) = [];
tGrid = tau0+1:T;

% --- 替换为浅灰、深灰、黑 的灰度配色 ---
CF = [
    0.7 0.7 0.7   % 1. 浅灰 (Small)
    0.4 0.4 0.4   % 2. 深灰 (Medium)
    0.0 0.0 0.0   % 3. 黑色 (Large)
];

scaleName = {'Small','Medium','Large'};
baseDir = fullfile(fileparts(thisDir), 'figures');

% 如果figures文件夹不存在，自动创建
if ~exist(baseDir, 'dir')
    mkdir(baseDir);
end

figPos = [80 80 250 200];     % 设置单张图的尺寸

xTickVals = [tau0 60 70 80 90 100];
xTickLabs = {'\tau_0','60','70','80','90','100'};

% ---------------- 1. beta estimation error ----------------
figure;
set(gcf, 'Position', figPos);
plot(tGrid, betaErr_FUCB(:,1), 'Color', CF(1,:), 'LineStyle', '-', 'LineWidth', 1.5); hold on;
plot(tGrid, betaErr_FUCB(:,2), 'Color', CF(2,:), 'LineStyle', '--', 'LineWidth', 1.5);
plot(tGrid, betaErr_FUCB(:,3), 'Color', CF(3,:), 'LineStyle', '-.', 'LineWidth', 1.5);
xlabel('Time index t');
xticks(xTickVals);
xticklabels(xTickLabs);
ylim([0,1.2])
lgd = legend(scaleName, 'Location', 'northeast'); 
lgd.ItemTokenSize = [15,8];
box on;

fullPath = fullfile(baseDir, 'one_beta_err.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);

% ---------------- 2. cumulative regret ----------------
figure;
set(gcf, 'Position', figPos);
plot(tGrid, regret_FUCB(:,1), 'Color', CF(1,:), 'LineStyle', '-', 'LineWidth', 1.5); hold on;
plot(tGrid, regret_FUCB(:,2), 'Color', CF(2,:), 'LineStyle', '--', 'LineWidth', 1.5);
plot(tGrid, regret_FUCB(:,3), 'Color', CF(3,:), 'LineStyle', '-.', 'LineWidth', 1.5);
xlabel('Time index t');
xticks(xTickVals);
xticklabels(xTickLabs);
box on;
lgd = legend(scaleName, 'Location', 'northwest'); 
lgd.ItemTokenSize = [15,8];
ylim([0,300])

fullPath = fullfile(baseDir, 'one_regret.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);

% ---------------- 3. test statistic ----------------
figure;
set(gcf, 'Position', figPos);
plot(tGrid, WPath_FUCB(:,1), 'Color', CF(1,:), 'LineStyle', '-', 'LineWidth', 1.5); hold on;
plot(tGrid, WPath_FUCB(:,2), 'Color', CF(2,:), 'LineStyle', '--', 'LineWidth', 1.5);
plot(tGrid, WPath_FUCB(:,3), 'Color', CF(3,:), 'LineStyle', '-.', 'LineWidth', 1.5);
xlabel('Time index t');
xticks(xTickVals);
xticklabels(xTickLabs);
box on;
lgd = legend(scaleName, 'Location', 'northwest'); 
lgd.ItemTokenSize = [15,8];

fullPath = fullfile(baseDir, 'one_test_stat.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);

% ---------------- 4. runtime (Boxplot) ----------------
figure;
set(gcf, 'Position', figPos);

% 浅灰、深灰、黑色
CF = [0.75 0.75 0.75;   % light gray
      0.35 0.35 0.35;   % dark gray
      0.00 0.00 0.00];  % black

boxplot(runtimeMat, ...
    'Labels', scaleName, ...
    'Colors', CF, ...
    'Symbol', 'o', ...
    'Widths', 0.5);

% 加粗所有 boxplot 线条
set(findobj(gca, 'type', 'line'), 'LineWidth', 1.2);

% 设置 outlier 为实心圆，并按 Small / Medium / Large 分别着色
hOut = findobj(gca, 'Tag', 'Outliers');

for i = 1:numel(hOut)
    x = get(hOut(i), 'XData');
    groupIdx = round(nanmedian(x));   % 对应第几个 box

    if groupIdx >= 1 && groupIdx <= size(CF, 1)
        set(hOut(i), ...
            'Marker', 'o', ...
            'MarkerFaceColor', CF(groupIdx, :), ...
            'MarkerEdgeColor', CF(groupIdx, :), ...
            'MarkerSize', 2, ...
            'LineStyle', 'none');
    end
end

box on;

fullPath = fullfile(baseDir, 'one_runtime.eps');
exportgraphics(gcf, fullPath, 'Resolution', 300);