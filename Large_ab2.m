clc;
clear;
close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(thisDir));

rng('default');

% ================================================================
% Part 1. User-specified model ingredients
% ================================================================
M = 500;
d = 50;
sigma2 = 0.01;
rho = 0.6;
Sigma = toeplitz(rho .^ (0:d-1));%
% Sigma = eye(d);
Phi = build_basis_matrix(M, d);
figure;plot(1:M,Phi);
tau0 = 50; % random set as T/5 ~ T/10
T = 1000;
w = 20; % w<tau0
initialRandomSteps = 20; % set as window size
m = 25; % set as sTrue*(2~3)

nRep = 500; % 500
nRepNull = 1000; % 1000
shiftMagnitudeSet = 0:0.2:1.8;
% parpool(20);

methodList = {'ORACLE','SB-FQCD','FUCB','SRAS','TS','CMAB'}; %SFUCB
nMethod = numel(methodList);
nShift = numel(shiftMagnitudeSet);

% ================================================================
% Part 2. Visual sanity check: normal mean vs abnormal mean
% ================================================================
a = 1;
I = [200 220];
shift = shiftMagnitudeSet(end);
[Yic, Yoc, betaStar, gTrue] = OCgen(Phi, Sigma, sigma2, shift, tau0, T,  'block', a, I);


mu0 = zeros(M, 1);
mu1 = shift * gTrue;

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

parpool(20)
calibrateThreshold = 1;
shift = 0;
targetARL0 = 200;
thresholdVec = nan(nMethod, 1);
ARL0 = zeros(nMethod, 1);
SDRL0 = zeros(nMethod, 1);
thresholdVec(:) = Inf;

if calibrateThreshold
    tic;
    Test_oracle = zeros(nRepNull,T);
    parfor s = 1:nRepNull
        [~,Test_oracle(s,:)] = sparse_bandit(Phi, Sigma, sigma2, shift, tau0, T, m, w, ...
            Inf,  'oracle', initialRandomSteps, 'block', a, I);
    end
    [thresholdVec(1),ARL0(1),SDRL0(1)] = FindCL(Test_oracle,targetARL0);
    toc;


    tic;
    Test_proposed = zeros(nRepNull,T);
    parfor s = 1:nRepNull
        [~,Test_proposed(s,:)] = sparse_bandit(Phi, Sigma, sigma2, shift, tau0, T, m, w, ...
            Inf,  'proposed', initialRandomSteps, 'block', a, I);
    end
    [thresholdVec(2),ARL0(2),SDRL0(2)] = FindCL(Test_proposed,targetARL0);
    toc;

    tic;
    Test_FUCB = zeros(nRepNull,T);
    parfor s = 1:nRepNull
        [~,Test_FUCB(s,:)] = sparse_bandit(Phi, Sigma, sigma2, shift, tau0, T, m, w, ...
            Inf,  'FUCB', initialRandomSteps, 'block', a, I);
    end
    [thresholdVec(3),ARL0(3),SDRL0(3)] = FindCL(Test_FUCB,targetARL0);
    toc;

    tic;
    Test_SRAS = zeros(nRepNull,T);
    parfor s = 1:nRepNull
        [Yic,Yoc] = OCgen(Phi, Sigma, sigma2, shift, tau0, T, 'block', a, I);
        [~,Test_SRAS(s,:)] = SRAS(Yic,Yoc,Inf,m,0);
    end
    [thresholdVec(4),ARL0(4),SDRL0(4)] = FindCL(Test_SRAS,targetARL0);
    toc;

    tic;
    Test_TS = zeros(nRepNull,T);
    parfor s = 1:nRepNull
        [Yic,Yoc] = OCgen(Phi, Sigma, sigma2, shift, tau0, T, 'block', a, I);
        [~,Test_TS(s,:)] = thom666(Yoc,Phi,sigma2,m,Inf,0);
    end
    [thresholdVec(5),ARL0(5),SDRL0(5)] = FindCL(Test_TS,targetARL0);
    toc;

    

    tic;
    Test_CMAB = zeros(nRepNull,T);
    parfor s = 1:nRepNull
        [Yic,Yoc] = OCgen(Phi, Sigma, sigma2, shift, tau0, T, 'block', a, I);
        [~,Test_CMAB(s,:)] = MABind(Yic,Yoc,m,Inf,0);
    end
    [thresholdVec(6),ARL0(6),SDRL0(6)] = FindCL(Test_CMAB,targetARL0);
    toc;

    save([thisDir,'\0ARL_save\Large_ab2\UCL'],'thresholdVec','ARL0','SDRL0');  
else
    load([thisDir,'\0ARL_save\Large_ab2\UCL']);
end


% ================================================================
% Oracle. Sweep shift magnitude and record detection delay
% ================================================================
for i = 2:nShift
    tic;
    shift = shiftMagnitudeSet(i);
    detectionDelayVec = nan(nRep, 1);
    Test_oracle_oc = nan(nRep, T);
    parfor s = 1:nRep
        [detectionDelayVec(s),Test_oracle_oc(s,:)] = sparse_bandit(Phi, Sigma, sigma2, shift, tau0, T, m, w, ...
            thresholdVec(1),  'oracle', initialRandomSteps, 'block', a, I);
    end

    meanDelay = mean(detectionDelayVec);
    stdDelay = std(detectionDelayVec);
    fname = [thisDir,'\0ARL_save\Large_ab2\ORACLE_','shift',num2str(i)];
    save(fname,'meanDelay','stdDelay');  
    toc;
end






% ================================================================
% Proposed. Sweep shift magnitude and record detection delay
% ================================================================
for i = 2:nShift
    tic;
    shift = shiftMagnitudeSet(i);
    detectionDelayVec = nan(nRep, 1);
    Test_proposed_oc = nan(nRep, T);
    parfor s = 1:nRep
        [detectionDelayVec(s),Test_proposed_oc(s,:)] = sparse_bandit(Phi, Sigma, sigma2, shift, tau0, T, m, w, ...
            thresholdVec(2),  'proposed', initialRandomSteps, 'block', a, I);
    end

    meanDelay = mean(detectionDelayVec);
    stdDelay = std(detectionDelayVec);
    fname = [thisDir,'\0ARL_save\Large_ab2\SB-FQCD_','shift',num2str(i)];
    save(fname,'meanDelay','stdDelay'); 
    toc;
end


% ================================================================
% FUCB without LASSO. Sweep shift magnitude and record detection delay
% ================================================================
for i = 2:nShift
    tic;
    shift = shiftMagnitudeSet(i);
    detectionDelayVec = nan(nRep, 1);
    Test_FUCB_oc = nan(nRep, T);
    parfor s = 1:nRep
        [detectionDelayVec(s),Test_FUCB_oc(s,:)] = sparse_bandit(Phi, Sigma, sigma2, shift, tau0, T, m, w, ...
            thresholdVec(3),  'FUCB', initialRandomSteps, 'block', a, I);
    end

    meanDelay = mean(detectionDelayVec);
    stdDelay = std(detectionDelayVec);
    fname = [thisDir,'\0ARL_save\Large_ab2\FUCB_','shift',num2str(i)];
    save(fname,'meanDelay','stdDelay');  
    toc;
end




% ================================================================
% SRAS. Sweep shift magnitude and record detection delay
% ================================================================
for i = 2:nShift
    tic;
    shift = shiftMagnitudeSet(i);
    detectionDelayVec = nan(nRep, 1);
    Test_SRAS_oc = nan(nRep, T);
    parfor s = 1:nRep
        [Yic, Yoc] = OCgen(Phi, Sigma, sigma2, shift, tau0, T, 'block', a, I);
        [detectionDelayVec(s),Test_SRAS_oc(s,:)] = SRAS(Yic,Yoc,thresholdVec(4),m,tau0);
    end

    meanDelay = mean(detectionDelayVec);
    stdDelay = std(detectionDelayVec);
    fname = [thisDir,'\0ARL_save\Large_ab2\SRAS_','shift',num2str(i)];
    save(fname,'meanDelay','stdDelay');  
    toc;
end

% ================================================================
% Tompson sampling with spike-slab. Sweep shift magnitude and record detection delay
% ================================================================
for i = 2:nShift
    tic;
    shift = shiftMagnitudeSet(i);
    detectionDelayVec = nan(nRep, 1);
    Test_TS_oc = nan(nRep, T);
    parfor s = 1:nRep
        [Yic, Yoc] = OCgen(Phi, Sigma, sigma2, shift, tau0, T, 'block', a, I);
        [detectionDelayVec(s),Test_TS_oc(s,:)] = thom666(Yoc,Phi,sigma2,m,thresholdVec(5),tau0);
    end

    meanDelay = mean(detectionDelayVec);
    stdDelay = std(detectionDelayVec);
    fname = [thisDir,'\0ARL_save\Large_ab2\TS_','shift',num2str(i)];
    save(fname,'meanDelay','stdDelay');  
    toc;
end



% ================================================================
% CMAB without LASSO. Sweep shift magnitude and record detection delay
% ================================================================
for i = 2:nShift
    tic;
    shift = shiftMagnitudeSet(i);
    detectionDelayVec = nan(nRep, 1);
    Test_CMAB_oc = nan(nRep, T);
    parfor s = 1:nRep
        [Yic, Yoc] = OCgen(Phi, Sigma, sigma2, shift, tau0, T, 'block', a, I);
        [detectionDelayVec(s),Test_CMAB_oc(s,:)] = MABind(Yic,Yoc,m,thresholdVec(6),tau0);
    end

    meanDelay = mean(detectionDelayVec);
    stdDelay = std(detectionDelayVec);
    fname = [thisDir,'\0ARL_save\Large_ab2\CMAB_','shift',num2str(i)];
    save(fname,'meanDelay','stdDelay');  
    toc;
end


delete(gcp('nocreate'))




