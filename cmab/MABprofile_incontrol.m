clear;
clc;
V=csvread('bases.csv');
d=size(V,2);
M=size(V,1);

m=20;
limit=52;
lambda_weight=0.9;
m00=0;
loss=@(t)sqrt(max(0,log((1-(1-lambda_weight).^(t))./lambda_weight)));
sigma_e2=0.05;
Lambda=eye(d);
T=1000;
T_nor=250;
T0=20;
delta=0;
nrep=500;
% offline
Y_nor=zeros(M,T_nor);
xi_nor=zeros(d,T_nor);
for t = 1:T_nor
xi_nor(:,t) = mvnrnd(zeros(d,1),Lambda,1);
Y_nor(:,t) = V*xi_nor(:,t)+ sqrt(sigma_e2)*randn(M,1);
end
normal_mean=mean(Y_nor,2);
normal_std=std(Y_nor,0,2);
covMat=1/(T_nor-1)*(Y_nor-repmat(normal_mean,1,T_nor))*(Y_nor-repmat(normal_mean,1,T_nor))';
    
Test = zeros(nrep,T);
Runlength = zeros(nrep,1);
parpool('local',8);
parfor ss=1:nrep

    %online
    Y=[];
    xi=zeros(d,T);
    mu=zeros(d,1);
    index = randsample(1:d,1);
    mu(index) = delta;
    for t = 1:T0
    xi(:,t) = mvnrnd(zeros(d,1),Lambda,1);
    Y(:,t) = V*xi(:,t)+ sqrt(sigma_e2)*randn(M,1);
    end
    for t = T0+1:T
    xi(:,t) = mvnrnd(mu,Lambda,1);
    Y(:,t) = V*xi(:,t)+ sqrt(sigma_e2)*randn(M,1);
    end

%     Y=(Y-repmat(normal_mean,1,T))./normal_std;

ss
[count1,Q,meanU1,runlen,test,index] = MABind(Y',m00,m,covMat,loss,lambda_weight,limit,T0);
% plot(test(20:T))
Test(ss,:)=test;
Runlength(ss) = runlen;
end
arl= mean(Runlength);
sdrl=std(Runlength);
delete(gcp('nocreate'))