clc;
clear;

n = 100; % Image size n*n
T = 1000; % The Total time length
bd = 4; % How smooth the defect is (Order of Spline)
% k0 = 20; % scale of the background
% bd0 = 3;
nq0=3;
k1  = 4; % scale of the defect
ndefect = 1;  % How many defect in the system
delta =0; % Defect magnitude
sigma = 0.05; % The noise level
Tin = 0; % When the change occur



repNum =1000;
tau=0;
m =30; % budget of pixels that can be observed for each image
lambda=0.1;
limit=1500;
m0=0;
% p=n*n;

% for i = 1:n
% for j = i+1:n
%     covMat(i,j) = 0.5^(abs(i-j));
%     covMat(j,i) = covMat(i,j);
% end
% end

loss = @(t) sqrt(max(0,log((1-(1-lambda).^t)./lambda)));
runlen=zeros(1,repNum);
Test = zeros(repNum,T-m0);
count1 = zeros(repNum,n);
meanU = count1;

parpool('local',8);
parfor rep =  1:repNum
    %generate the background
    B0 = fourierbasis(nq0,n);
    B00 = B0;
    mu0 = zeros(1,size(B00,2));
    sigma0 = 0.1;
    Sigma0 = eye(size(B00,2),size(B00,2))*sigma0;
    theta0=zeros(T,size(B00,2));
    % sample IC data
    Y = [];
%     for i = 1:T
%     theta0(i,:) = mvnrnd(mu0,Sigma0,1);
%     Y(:,:,i) = reshape(B00*theta0(i,:)',n,n)+ sigma*randn(n);
%     end
    for i = 1:T
    theta0(i,:) = mvnrnd(mu0,Sigma0,1);
    Y(:,i) = B00*theta0(i,:)'+ sigma*randn(n,1);
    end
    normal_mean=mean(Y,2);
    normal_std=std(Y,0,2);
    covMat = eye(n)*sigma^2+sigma0.*B00*eye(size(B00,2))*B00'; 

    
    %Generate Defect by Spline Bases (We can change this later if we need to generate other shape of defect
    knots = [ones(1,bd) linspace(1,n,round(n/k1)) n * ones(1,bd)];
    nKnots = length(knots) - (bd-1);
    kspline = spmak(knots,eye(nKnots));
    kd = fnder(kspline);
    B=spval(kspline,1:n)';
    B = B(:,3:end-2);
    
    
    kcoef = size(B,2);
    Z = zeros(kcoef,1);

    for i = 1:ndefect
        Z(randsample(1:kcoef,1)) = 1;
    end
    %imagesc(Z)

    D = B*Z;
    %imagesc(D)
    %colorbar;
    for i = (Tin+1):T
        Y(:,i) = Y(:,i)+delta*D;
    end
    Y=(Y-repmat(normal_mean,1,T))./normal_std;

    rep
   %tic;
    [count1,Q,meanU,runlen(rep),Test(rep,:)] = MABind(Y',m0,m,covMat,loss,lambda,limit,tau);
%     [count1,Q,meanU,runlen,Test] = MABind(Y',m0,m,covMat,loss,lambda,limit,tau);
%     plot(Test)
    %toc;
end
arl_mab= mean(runlen);
sdrl_mab= std(runlen);
delete(gcp('nocreate'))
% limit = 2370;
% for rep = 1:repNum
%     ind =  find(Test(rep,tau+1:end)>limit,1);
%     if isempty(ind)
%         runlen1(rep) = T;
%     else
%         runlen1(rep) = ind;
%     end
% end
% mean(runlen1)
% std(runlen1)
% runlen(find(runlen<=tau)) = NaN;
% nanmean(runlen,2)
% 
% limit = 12.15
% for rep = 1:repNum
%     ind =  find(Test(rep,:)>limit,1);
%     if isempty(ind)
%         runlen1(rep) = k0-m0;
%     else
%         runlen1(rep) = ind;
%     end
% end
% mean(runlen1)
% std(runlen1)
% 
% for i = 1:p
% figure('position',[40,100,1500,1000]);
%     plot(WU(:,i),'r');
% hold on;
% plot(WL(:,i),'k');
% plot(W0(:,i));
% end
% 
% % power
% ll = distparams1{1}*inv(covMat)*distparams1{1}';
% ncx = ncx2cdf(limit,p,ll);
% OCARL = 1/(1-ncx);