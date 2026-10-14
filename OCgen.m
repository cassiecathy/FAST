function [Yic, Yoc, betaStar, gTrue] = OCgen(Phi, Sigma, sigma2, shift, tau0, T, mode, varargin)

[M, d] = size(Phi);

betaStar = zeros(d,1);
gTrue = zeros(M,1);

switch lower(mode)
    case 'beta'
        betaSupport = varargin{1};
        betaSigns = varargin{2};

        for k = 1:numel(betaSupport)
            betaStar(betaSupport(k)) = shift * betaSigns(k);
        end

    case 'block'
        a = varargin{1};
        I = varargin{2};

        if numel(I) == 2 && all(mod(I,1) == 0)
            idx = I(1):I(2);
        else
            idx = I(:)';
        end

        gTrue(idx) = a;

    otherwise
        error('mode must be ''beta'' or ''block''.')
end

Yic = zeros(M,T);
Yoc = zeros(M,T);

for t = 1:T
    x_ic = mvnrnd(zeros(d,1), Sigma).';
    Yic(:,t) = Phi * x_ic + sqrt(sigma2) * randn(M,1);
end

for t = 1:tau0
    x_oc = mvnrnd(zeros(d,1), Sigma).';
    Yoc(:,t) = Phi * x_oc + sqrt(sigma2) * randn(M,1);
end

for t = tau0+1:T
    x_oc = mvnrnd(zeros(d,1), Sigma).';

    if strcmpi(mode,'beta')
        Yoc(:,t) = Phi * (x_oc + betaStar) + sqrt(sigma2) * randn(M,1);
    else
        Yoc(:,t) = Phi * x_oc + shift * gTrue + sqrt(sigma2) * randn(M,1);
    end
end

if strcmpi(mode,'block')

    muPre = mean(Yoc(:,1:tau0),2); % M*1
    muPost = mean(Yoc(:,tau0+1:end),2); % M*1
    deltaHat = muPost - muPre;
    betaStar = Phi \ deltaHat;
end

end