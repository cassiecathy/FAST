clc;
clear;


V=csvread('bases.csv');

d=size(V,2);
p=size(V,1);

sigma_e2=0.03;
Lambda=eye(d);
T_nor=2000;
T=1000;
T0=20;
m=20; % budget of pixels that can be observed for each image
lambda=0.1;
h=1.1;
delta1=0.1;
delta=0;
cl=44;
nrep=100;

X_nor=zeros(p,T_nor);
xi_nor=zeros(d,T_nor);
for t = 1:T_nor
xi_nor(:,t) = mvnrnd(zeros(d,1),Lambda,1);
X_nor(:,t) = V*xi_nor(:,t)+ sqrt(sigma_e2)*randn(p,1);
end
normal_mean=mean(X_nor,2);
normal_std=std(X_nor,0,2);
X_nor=X_nor-normal_mean;
    

Test = zeros(nrep,T);
Runlength = zeros(nrep,1);
parpool('local',64);
parfor ss=1:nrep

    %online
    X=[];
    xi=zeros(d,T);
    mu=zeros(d,1);
    index = randsample(1:d,1);
    mu(index) = delta;
    for t = 1:T0
        xi(:,t) = mvnrnd(zeros(d,1),Lambda,1);
        X(:,t) = V*xi(:,t)+ sqrt(sigma_e2)*randn(p,1);
    end
    for t = T0+1:T
        xi(:,t) = mvnrnd(mu,Lambda,1);
        X(:,t) = V*xi(:,t)+ sqrt(sigma_e2)*randn(p,1);
    end
     X=X-normal_mean;
    
    ss
    [Q,runlen,index8] = SRAS(X,X_nor,h,delta1,cl,lambda,m,T0);
% plot(Q(21:1000))
Test(ss,:)=Q;
Runlength(ss) = runlen;
end
arl= mean(Runlength);
sdrl=std(Runlength);
delete(gcp('nocreate'))

% 
% for t=21:50
%     plot(Y(:,t));
%     pause(0.8)
% end



