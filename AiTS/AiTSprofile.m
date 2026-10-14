clc;
clear;


V=csvread('bases.csv');

d=size(V,2);
p=size(V,1);

sigma_e2=0.05;
Lambda=eye(d);
T_nor=250;
T=1000;
T0=100;
m=20; % budget of pixels that can be observed for each image
m0=0;
g=repmat(1/p,p,1);
k=1.5;
h=0.825;


nrep=200;

Y_nor=zeros(p,T_nor);
xi_nor=zeros(d,T_nor);
for t = 1:T_nor
xi_nor(:,t) = mvnrnd(zeros(d,1),Lambda,1);
Y_nor(:,t) = V*xi_nor(:,t)+ sqrt(sigma_e2)*randn(p,1);
end
normal_mean=mean(Y_nor,2);
normal_std=std(Y_nor,0,2);
    

Test = zeros(nrep,T);
Runlength = zeros(nrep,1);

delta1set=0.6:0.2:2.4;
arl=zeros(size(delta1set,2),1);
stdrl=zeros(size(delta1set,2),1);

parpool('local',8)
for dd=1:size(delta1set,2)
     delta=delta1set(dd);
parfor ss=1:nrep

    %online
    Y=[];
    xi=zeros(d,T);
    mu=zeros(d,1);
    index = randsample(1:d,1);
    mu(index) = delta;
    for t = 1:T0
        xi(:,t) = mvnrnd(zeros(d,1),Lambda,1);
        Y(:,t) = V*xi(:,t)+ sqrt(sigma_e2)*randn(p,1);
    end
    for t = T0+1:T
        xi(:,t) = mvnrnd(mu,Lambda,1);
        Y(:,t) = V*xi(:,t)+ sqrt(sigma_e2)*randn(p,1);
    end
    Y1=(Y-normal_mean)./normal_std;
   dd 
    ss
    [test,runlen,index8] = AiTS(Y1',g,m,m0,k,h,T0);
% plot(test(101:1000))
Test(ss,:)=test;
Runlength(ss) = runlen;
end
arl(dd)= mean(Runlength);
stdrl(dd)= std(Runlength);
end
delete(gcp('nocreate'))

% 
% for t=21:50
%     plot(Y(:,t));
%     pause(0.8)
% end



