clc;
clear;


V=csvread('bases.csv');


d=size(V,2);
M=size(V,1);

figure; plot(1:M, V, 'LineWidth', 1);

% theta=zeros(kb,T_nor);
% e=zeros(p,T_nor);
% for t=1:T_nor
%     theta(:,t)=inv(B00'*B00)*B00'*Y0(:,t);
%     e(:,t)=(eye(p)-B00*inv(B00'*B00)*B00')*Y0(:,t);
% end
% Sigma=cov(theta');
% sigma00=mean(diag(Sigma));
% sigma=sqrt(sum((e-mean(e,2)).^2,'all')/(p*T_nor));

sigma00=1;
sigma_e2=0.05;
Lambda=eye(d);
T=1000;
T0=20;
m=20; % budget of pixels that can be observed for each image
niter=50;
controlimit=1.44*10^-3;
lambda0=0.9;
B=V;
w=0.2*ones(d,1);
sigma9=2*ones(d,1);

nrep=1000;

deltaset=0.6:0.2:2.4;

arl=zeros(size(deltaset,2),1);
sdrl=zeros(size(deltaset,2),1);

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
    dd
    ss
[~,index8,theta8,theta80,test,runlen] = thom666(Y',V,B,m,sqrt(sigma_e2),niter,w,sigma9,lambda0,controlimit,sigma00,T0);

Test(ss,:)=test;
Runlength(ss) = runlen;
end
arl(dd)= mean(Runlength);
sdrl(dd)=std(Runlength);
end
delete(gcp('nocreate'))
