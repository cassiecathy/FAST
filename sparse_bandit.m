function [detectionDelay, WPath, betaErrPath, regretPath, ZHistory, timePerStepPath] = ...
    sparse_bandit(Phi, Sigma, sigma2, shift, tau0, T, m, w, threshold, ...
    method, initialRandomSteps, mode, varargin)
%SIMULATE_SPARSE_BANDIT_ONCE One online replication under a chosen sensing method.
% Detection delay is reported only when tauHat >= tau0. Early alarms are treated
% as false alarms and excluded from delay averaging via detectionDelay = NaN.

% Dimensions
[M, d] = size(Phi);
[~, yFull, betaStar] = OCgen(Phi, Sigma, sigma2, shift, tau0, T, mode, varargin{:});

if shift == 0
    tau0 = 0;
    betaStar = zeros(d,1);
end

sigma = sqrt(sigma2);
sTrue = nnz(abs(betaStar) > 1e-12);

% Fixed numerical constants
maxCdIter = 100;
cdTol = 1e-6;
delta = 0.05;
ucbBetaTol = 1e-10;
maxOracleComb = 50000;
phiFloor = 0.05;

Omega = Phi * Sigma * Phi.' + sigma2 * eye(M);
Omega = (Omega + Omega.') / 2;

[Istar, Zstar] = oracle_subset_information(betaStar, Phi, Omega, m, maxOracleComb);

S = zeros(d, d);
b = zeros(d, 1);
AHistory = cell(T, 1);
gHistory = cell(T, 1);
betaPrev = zeros(d, 1);
W = 0;
cumRegret = 0;

WPath = nan(1, T);
betaHatPath = nan(d, T);
betaErrPath = nan(1, T);
regretPath = nan(1, T);
overlapPath = nan(1, T);
timePerStepPath = nan(1, T);
ZHistory = cell(T, 1);


tauHat = T;

for t = 1:T
    stepClock = tic;

    % 1) Online beta update from data observed up to t-1.
    nWindow = min(t - 1, w);
    if nWindow == 0
        betaHat = zeros(d, 1);
    else
        if strcmpi(method, 'FUCB')
            betaHat = (S + 1e-8 * eye(d)) \ b;
    
        elseif any(strcmpi(method, {'proposed', 'random', 'oracle'}))
            lambda = min(4 * sigma * sqrt(2 * log(2 * d / delta) / nWindow),0.1);
            betaHat = lasso_from_sufficient_stats(S, b, nWindow, lambda, betaPrev, ...
                maxCdIter, cdTol);
        else
            error('Unknown method: %s', method);
        end
    end
    betaHatPath(:, t) = betaHat;
    if t <=tau0
        betaErrPath(t) = norm(betaHat, 2);
    else
        betaErrPath(t) = norm(betaHat - betaStar, 2) / max(norm(betaStar, 2), 1e-8);
    end
    

    % 2) Choose sensing subset.
    switch lower(method)
        case 'random'
            Zt = sort(randperm(M, m));
            Szt = subset_information(Zt, Phi, Omega);

        case 'oracle'
            Zt = Zstar;
            Szt = subset_information(Zt, Phi, Omega);

        otherwise % proposed
            if (norm(betaHat, 2) <= ucbBetaTol) || (t <= initialRandomSteps)   % if in-control or if very early
                Zt = sort(randperm(M, m));
                Szt = subset_information(Zt, Phi, Omega);
            else
                % Real-time phi estimate and theorem-driven alpha.
                Gphi = S / nWindow;
                phiConst = estimate_phi_const(Gphi, betaStar, betaHat, method, phiFloor);
                alpha0 = sqrt((72 * sigma2 / (phiConst^2)) * max(sTrue, 1) * log(2 * d / delta));
                alpha0 = 20*alpha0*(1/(1+exp(-50*(betaErrPath(t)-0.5))));
                Vt = safe_inverse(S);
                [Zt, Szt] = greedy_select_subset_ucb(betaHat, Vt, Phi, Omega, m, alpha0);

                % S_history{t} = S;
                % alpha0_history(t) = alpha0;
                % Vt_history{t} = Vt;
                % debug_show{1,t} = min(eig((S_history{t}+S_history{t}')/2));
                % debug_show{2,t} = norm(Vt_history{t}, 2);

            end
    end
    ZHistory{t} = Zt;
    overlapPath(t) = numel(intersect(Zt, Zstar)) / m;

    % 3) Generate current sample.
    yZ = yFull(Zt,t);

    reward = 0.5 * betaStar.' * Szt * betaStar;
    if t > tau0
        cumRegret = cumRegret + max(0, real(Istar - reward));
    end
    regretPath(t) = cumRegret;

    % 4) Detection statistic.
    ell = llr_increment(yZ, Zt, betaHat, Phi, Omega);
    W = max(0, W + ell);
    WPath(t) = W;

    
    if W >= threshold && t>tau0
        tauHat = t;
        break;
    end

    % 5) Update sliding-window sufficient statistics with current sample.
    [At, gt] = sufficient_stats_from_obs(yZ, Zt, Phi, Omega);
    AHistory{t} = At;
    gHistory{t} = gt;
    S = S + At;
    b = b + gt;

    tOld = t - w;
    if tOld >= 1
        S = S - AHistory{tOld};
        b = b - gHistory{tOld};
    end
    
    timePerStepPath(t) = toc(stepClock);
    betaPrev = betaHat;
end

detectionDelay = tauHat - tau0;


% tFinal = min(tauHat, T);
% WPath = WPath(1:tFinal);
% betaErrPath = betaErrPath(1:tFinal);
% regretPath = regretPath(1:tFinal);
% ZHistory = ZHistory(1:tFinal);
% betaHatPath = betaHatPath(:, 1:tFinal);
% overlapPath = overlapPath(1:tFinal);
% timePerStepPath = timePerStepPath(1:tFinal);


end

% ================================================================
% Local functions
% ================================================================
function phiConst = estimate_phi_const(Gphi, betaStar, betaHat, method, phiFloor)
    Gphi = (Gphi + Gphi.') / 2;

    switch lower(method)
        case 'oracle'
            Sidx = find(abs(betaStar) > 1e-12);
        otherwise
            Sidx = find(abs(betaHat) > 1e-6);
            if isempty(Sidx)
                Sidx = find(abs(betaStar) > 1e-12);
            end
    end

    if isempty(Sidx)
        phiConst = 1.0;
        return;
    end

    Gs = Gphi(Sidx, Sidx);
    eigMin = min(eig((Gs + Gs.') / 2));
    phiConst = sqrt(max(eigMin, phiFloor^2));
end

function betaHat = lasso_from_sufficient_stats(S, b, nWindow, lambda, betaStart, ...
    maxCdIter, cdTol)
    d = length(betaStart);
    betaHat = betaStart(:);
    if nWindow <= 0
        betaHat = zeros(d, 1);
        return;
    end

    for iter = 1:maxCdIter
        betaOld = betaHat;
        for j = 1:d
            partialResidual = (b(j) - S(j, :) * betaHat + S(j, j) * betaHat(j)) / nWindow;
            denom = S(j, j) / nWindow;
            betaHat(j) = soft_threshold(partialResidual, lambda) / max(denom, 1e-12);
        end
        if norm(betaHat - betaOld, inf) < cdTol
            break;
        end
    end
end

function value = soft_threshold(x, tau)
    value = sign(x) .* max(abs(x) - tau, 0);
end

function [Zbest, Sbest] = greedy_select_subset_ucb(betaHat, Vt, Phi, Omega, m, alpha0)
    M = size(Phi, 1);
    Zbest = [];
    Sbest = zeros(size(Phi, 2));
    available = 1:M;

    for k = 1:m
        bestScore = -inf;
        bestSensor = available(1);
        bestS = [];
        for idx = 1:length(available)
            candidateSet = sort([Zbest, available(idx)]);
            Scand = subset_information(candidateSet, Phi, Omega);
            exploit = 0.5 * betaHat.' * Scand * betaHat;
            gradReward = Scand * betaHat;
            bonus = alpha0 * sqrt(max(real(gradReward.' * Vt * gradReward), 0));
            score = real(exploit + bonus);
            if score > bestScore
                bestScore = score;
                bestSensor = available(idx);
                bestS = Scand;
            end
        end
        Zbest = sort([Zbest, bestSensor]);
        Sbest = bestS;
        available(available == bestSensor) = [];
    end
end

function [At, gt] = sufficient_stats_from_obs(yZ, Zt, Phi, Omega)
    PhiZ = Phi(Zt, :);
    invOmegaZ = safe_inverse(Omega(Zt, Zt));
    At = PhiZ.' * invOmegaZ * PhiZ;
    gt = PhiZ.' * invOmegaZ * yZ;
end

function ell = llr_increment(yZ, Zt, betaHat, Phi, Omega)
    PhiZ = Phi(Zt, :);
    invOmegaZ = safe_inverse(Omega(Zt, Zt));
    muZ = PhiZ * betaHat;
    ell = real(muZ.' * invOmegaZ * yZ - 0.5 * muZ.' * invOmegaZ * muZ);
end

function Ssubset = subset_information(Zt, Phi, Omega)
    PhiZ = Phi(Zt, :);
    invOmegaZ = safe_inverse(Omega(Zt, Zt));
    Ssubset = PhiZ.' * invOmegaZ * PhiZ;
end

function [Istar, Zstar] = oracle_subset_information(betaStar, Phi, Omega, m, maxOracleComb)
    M = size(Phi, 1);
    logNumComb = gammaln(M + 1) - gammaln(m + 1) - gammaln(M - m + 1);
    numComb = exp(logNumComb);

    if numComb <= maxOracleComb
        combos = nchoosek(1:M, m);
        bestValue = -inf;
        Zstar = combos(1, :);
        for rr = 1:size(combos, 1)
            Scand = subset_information(combos(rr, :), Phi, Omega);
            value = real(0.5 * betaStar.' * Scand * betaStar);
            if value > bestValue
                bestValue = value;
                Zstar = combos(rr, :);
            end
        end
        Istar = bestValue;
    else
        Vdummy = zeros(size(Phi, 2));
        [Zstar, Sstar] = greedy_select_subset_ucb(betaStar, Vdummy, Phi, Omega, m, 0);
        Istar = real(0.5 * betaStar.' * Sstar * betaStar);
    end
end

function Ainv = safe_inverse(A, tau)
    if nargin < 2
        tau = 1e-2;
    end
    A = (A + A.') / 2;
    [Q, D] = eig(A);
    lam = real(diag(D));
    lam = max(lam, tau);
    Ainv = Q * diag(1 ./ lam) * Q';
    Ainv = (Ainv + Ainv.') / 2;
end
