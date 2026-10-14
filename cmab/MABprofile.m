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
deltaset=1:0.2:1.4;
nrep=500;
arl=zeros(size(deltaset,2),1);
sdrl=zeros(size(deltaset,2),1);

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
    
parpool('local',8);
for dd=1:size(deltaset,2)
    delta=deltaset(dd);
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
dd
ss
[count1,Q,meanU1,runlen,test,index] = MABind(Y',m00,m,covMat,loss,lambda_weight,limit,T0);
Test(ss,:)=test;
Runlength(ss) = runlen;
end
arl(dd)= mean(Runlength);
sdrl(dd)=std(Runlength);
end
delete(gcp('nocreate'))