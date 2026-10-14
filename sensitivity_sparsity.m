clc;
clear;
close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(thisDir));

rng('default');

% ================================================================
% Part 1. User-specified model ingredients
% ================================================================
M = 100;
d = 20;
sigma2 = 0.01;
Sigma = eye(d);%toeplitz(rho .^ (0:d-1));
Phi = build_basis_matrix(M, d);
% figure;plot(1:M,Phi);
tau0 = 50; % random set as T/5 ~ T/10
T = 1000;
w = 15; % w<tau0
initialRandomSteps = 20; % set as window size
m = 10; % set as sTrue*(2~3)

nRep = 1000; % 500
nRepNull = 2000; % 1000
% parpool(20);

sparseList = [1,2,3,4,5]; %SFUCB
nSparse = numel(sparseList);
supportCell = {
    [10], ...
    [6 14], ...
    [4 10 16], ...
    [3 7 13 17], ...
    [2 6 10 14 18]
};

% ================================================================
% Part 2. Visual sanity check: normal mean vs abnormal mean
% ================================================================
betaSupport = [2 6 10 14 18];
betaSigns = [1 1 1 1 1]; % in OC_gen also
shift = 1;
[Yic, Yoc, betaStar] = OCgen(Phi, Sigma, sigma2, shift, tau0, T, 'beta', betaSupport, betaSigns);

mu0 = zeros(M, 1);
mu1 = Phi * betaStar;

figure;
plot(1:M, mu0, 'LineWidth', 1.8); hold on;
plot(1:M, mu1, 'LineWidth', 1.8);
legend('Normal mean', 'Abnormal mean', 'Location', 'best');
xlabel('Sensor index');
ylabel('Mean signal');
title('Normal and abnormal mean functions');
grid on;


% A few sample curves under H0 and H1
figure;
subplot(2,1,1); hold on; plot(1:M, Yic(:,1:3), 'LineWidth', 1.0);
title('Sample normal curves'); xlabel('Sensor index'); ylabel('Signal');
subplot(2,1,2); hold on; plot(1:M, Yoc(:,tau0:tau0+3), 'LineWidth', 1.0);
title('Sample abnormal curves'); xlabel('Sensor index'); ylabel('Signal');



% ================================================================
% Part 3. Threshold calibration for each method separately
% ================================================================
close all;

calibrateThreshold = 1;
shift = 0;
targetARL0 = 200;
thresholdVec = nan(nSparse, 1);
ARL0 = zeros(nSparse, 1);
SDRL0 = zeros(nSparse, 1);
thresholdVec(:) = Inf;

if calibrateThreshold
    for j = 1:nSparse
        tic;
        Test_proposed = zeros(nRepNull,T);
        sTrue = sparseList(j);
        betaSupport = supportCell{j};
        betaSigns = ones(1, sTrue);
        parfor s = 1:nRepNull
            [~,Test_proposed(s,:)] = sparse_bandit(Phi, Sigma, sigma2, shift, tau0, T, m, w, ...
                Inf,  'proposed', initialRandomSteps, 'beta', betaSupport, betaSigns);
        end
        [thresholdVec(j),ARL0(j),SDRL0(j)] = FindCL(Test_proposed,targetARL0);
        toc;
    end


    save([thisDir,'\0ARL_save\Sensitivity_sparsity_s\UCL'],'thresholdVec','ARL0','SDRL0');  
else
    load([thisDir,'\0ARL_save\Sensitivity_sparsity_s\UCL']);
end




% ================================================================
% Proposed. Sweep shift magnitude and record detection delay
% ================================================================
shiftList = 1 ./ sqrt(sparseList);
for j = 1:nSparse
    tic;
    detectionDelayVec = nan(nRep, 1);
    Test_proposed_oc = nan(nRep, T);
    sTrue = sparseList(j);
    betaSupport = supportCell{j};
    betaSigns = ones(1, sTrue);
    shift = 1 / sqrt(sTrue);
    parfor s = 1:nRep
        [detectionDelayVec(s),Test_proposed_oc(s,:)] = sparse_bandit(Phi, Sigma, sigma2, shift, tau0, T, m, w, ...
            thresholdVec(j),  'proposed', initialRandomSteps, 'beta', betaSupport, betaSigns);
    end

    meanDelay = mean(detectionDelayVec);
    stdDelay = std(detectionDelayVec);
    fname = [thisDir,'\0ARL_save\Sensitivity_sparsity_s\SB-FQCD_','sparse',num2str(j)];
    save(fname,'meanDelay','stdDelay'); 
    toc;
end



delete(gcp('nocreate'))




