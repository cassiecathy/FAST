clc;clear;

n=100;
censdim=30;% budget of pixels that can be observed for each image

shift = 0;
g=repmat(1/(n+1),n+1,1);
g(n+1)=0.1;
g(1:n)=repmat((1-g(n+1))/n,n,1);
delta=0.01;
k=1;
lambda = (n*mean(g(1:n))-(n-censdim)*delta)/(1-g(n+1)^(censdim/n));
lambda0 = g(n+1)^(1-censdim/n);
% kmax=min(lambda^2/g(1)-2*lambda,lambda0^2/g(n+1)-2*lambda0)+(n-censdim)*(delta^2/g(1)-2*delta)+1;

m0=0;
T = 1000; % The Total time length
bd = 4; % How smooth the defect is (Order of Spline)
% k0 = 20; % scale of the background
% bd0 = 3;
nq0 = 3;
k1  = 4; % scale of the defect
ndefect = 1;  % How many defect in the system
delta1 =0; % Defect magnitude
sigma = 0.05; % The noise level
Tin =0; % When the change occur
%cutoff = 0.5;

repNum =1000;
climit=30;
Tval = zeros(repNum, T);
id=zeros(repNum,1);
runlen=zeros(1,repNum);
delta1set=[0.5,0.8,1.0,1.5,2.0];
arl=zeros(1,size(delta1set,2));
stdrl=zeros(1,size(delta1set,2));

parpool('local',8)
for dd=1:size(delta1set,2)
     delta1=delta1set(dd);
  parfor tt = 1 : repNum
      %generate the background
%     knots0 = [ones(1,bd) linspace(1,n,round(n/k0)) n * ones(1,bd0)];
%     nKnots0 = length(knots0) - (bd0-1);
%     kspline0 = spmak(knots0,eye(nKnots0));
%     kd = fnder(kspline0);
%     B0=spval(kspline0,1:n)';
     B0 = fourierbasis(nq0,n);
    B00 = B0;
    mu0 = zeros(1,size(B00,2));
    sigma0 = 0.1;
    Sigma0 = eye(size(B00,2),size(B00,2))*sigma0;
    theta0 = zeros(T,size(B00,2));
    % sample IC data
    Y = [];
    for i = 1:T
    theta0(i,:) = mvnrnd(mu0,Sigma0,1);
    Y(:,i) = B00*theta0(i,:)'+ sigma*randn(n,1);
    end
    normal_mean=mean(Y,2);
    normal_std=std(Y,0,2);
    
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
    for i = (Tin+1):T
        Y(:,i) = Y(:,i)+delta1*D;
    end
    Y=(Y-repmat(normal_mean,1,T))./normal_std;
    
     dd
     tt
    [Tval(tt,:),runlen(tt)] = naCUSUM(-Y',censdim,m0,delta,k,climit,lambda,lambda0,g);   
  end
   arl(dd)= mean(runlen);
   stdrl(dd)= std(runlen);
end
 delete(gcp('nocreate'))
%% get control limit
% limit = 2.39
% for rep = 1:repNum
%     ind =  find(Tval(rep,2:end)>limit,1);
%     if isempty(ind)
%         runlen1(rep) = T;
%     else
%         runlen1(rep) = ind;
%     end
% end
% mean(runlen1)
% std(runlen1)

