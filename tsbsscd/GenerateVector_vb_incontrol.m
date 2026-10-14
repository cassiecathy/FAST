clc;
clear;

n = 100; % Image size n*n
T = 50; % The Total time length
bd = 5; % How smooth the defect is (Order of Spline)
k0 = 20; % scale of the background
bd0 = 3;
nq0 = 3;
k1  = 2; % scale of the defect
ndefect = 1;  % How many defect in the system
delta =1; % Defect magnitude
sigma = 0.05; % The noise level
Tin = 1; % When the change occur


m=100; % budget of pixels that can be observed for each image
niter=50;
controlimit=2.2*10^-2;
controlimit=inf;
lambda0=0.9;
nrep_oc=1000;


parpool('local',8);
parfor ss=1:nrep_oc
    %generate the background
    knots0 = [ones(1,bd) linspace(1,n,round(n/k0)) n * ones(1,bd0)];
    nKnots0 = length(knots0) - (bd0-1);
    kspline0 = spmak(knots0,eye(nKnots0));
    kd = fnder(kspline0);
    B0=spval(kspline0,1:n)';
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
        Y(:,i) = Y(:,i)+delta*D;
    end
     Y1 = Y;
    nq = size(B,2); 
    w=0.9*ones(nq,1);
    sigma9=3*ones(nq,1);
    sigma00=0.3;
    ss
    warning off;

[~,index8,theta8,theta80,Test,runlength,time] = thom666(Y1',B00,B,m,sigma,niter,w,sigma9,lambda0,controlimit,sigma00);
%  [~,index8,theta8,theta80,Test,runlength] = thom7(Y1',B00,B,m,sigma,controlimit,sigma0,u,a,s2);
TTest(ss,:)=Test;
RRunlength(ss) = runlength;
end
ARL_oc= mean(RRunlength);
sdrl_oc=std(RRunlength);
delete(gcp('nocreate'))
