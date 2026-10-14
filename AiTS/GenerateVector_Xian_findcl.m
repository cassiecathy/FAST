clc;
clear;
n = 30; % Image size n*n
T = 2000; % The Total time length
bd = 5; % How smooth the defect is (Order of Spline)
k0 = 20; % scale of the background
bd0 = 3;
nq0 = 2;
k1  = 2; % scale of the defect
ndefect = 1;  % How many defect in the system
delta1 =0; % Defect magnitude
sigma = 0.05; % The noise level
Tin = 1; % When the change occur
Y = zeros(n,n,T); % Generated Videos
%cutoff = 0.5;
g = repmat(1/(n+1),n+1,1);
delta = 0.07;
k = 1;
censdim=n;% budget of pixels that can be observed for each image
lambda = (n*mean(g(1:n))-(n-censdim)*delta*0.53)/(1-g(n+1)^(censdim/n));
lambda0 = g(n+1)^(1-censdim/n)*3.5;

repNum=2000;
climit= inf;
r= censdim;
m0=0;
Tval = zeros(repNum, T);
id=zeros(repNum,1);
runlen=zeros(1,repNum);
ARL0=200;
parpool('local',12);
parfor ss=1:repNum
    knots0 = [ones(1,bd) linspace(1,n,round(n/k0)) n * ones(1,bd0)];
    nKnots0 = length(knots0) - (bd0-1);
    kspline0 = spmak(knots0,eye(nKnots0));
    kd = fnder(kspline0);
    B0=spval(kspline0,1:n)';
%    B0 = B0(:,4:end-2);
 %   B0 = B0(:,5:end-2);
     B0 = fourierbasis(nq0,n);
    B00 = B0;
   % B00 = kron(B0,B0);
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
    %colorbar;

    for i = (Tin+1):T
        Y(:,i) = Y(:,i)+delta1*D;
    end
 %   Y1 = reshape(Y,n*n,T);
     Y1 = Y;
    ss
    [Tval(ss,:),runlen(ss)] = naCUSUM(Y1',censdim,m0,delta,k,climit,lambda,lambda0,g); 
% TTest=Test;
% RRunlength=runlength;
TTest(ss,:) = Tval(ss,:);
RRunlength(ss) = runlen(ss);
end

% TTest=[TTest3;TTest4];
% RRunlength=[RRunlength3 RRunlength4];
%TTest(:,1)=[];
%TTest=TTest*1000;
CLmax=max(TTest(:));
CLmin=0;
cl=CLmax;
bia=1;
ARL(1)=0;
ARL(2)=mean(RRunlength);
TTest1=TTest';
i=2;
while(abs(ARL(i)-ARL0)>bia && abs(CLmax-CLmin)>10^-5)
    if ARL(i)>ARL0
        CLmax=cl;
        cl=0.5*(CLmax+CLmin);
    else
        CLmin=cl;
        cl=0.5*(CLmax+CLmin);
    end
    i=i+1;
    TTest2=TTest1;
    TTest2(abs(TTest1)<=cl)=0;
    [u,v]=find(TTest2);
    ARL(i)=mean(u(diff([0;v])~=0));
    ARL(i)
end
fprintf('ARL of in-control process is %f\n',ARL(i))
fprintf('The controlimit is %f',cl)
save('GenerateVector_Xian_findcl1.mat')
delete(gcp('nocreate'))