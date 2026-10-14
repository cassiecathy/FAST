clc;
clear;

n = 100; % Image size n*n
T = 1000; % The Total time length
bd = 5; % How smooth the defect is (Order of Spline)
k0 = 20; % scale of the background
bd0 = 3;
nq0 = 3;
k1  = 2; % scale of the defect
ndefect = 1;  % How many defect in the system
delta =0; % Defect magnitude
sigma = 0.05; % The noise level
Tin = 1; % When the change occur


m=5; % budget of pixels that can be observed for each image
niter=50;
controlimit=2.2*10^-2;
lambda0=0.9;
nrep_oc =1000;
deltaset=0.1:0.1:1;
ARL_oc=zeros(1,size(deltaset,2));
sdrl_oc=zeros(1,size(deltaset,2));
%E=zeros(1,T);

parpool('local',8);
for dd=1:size(deltaset,2)
    delta=deltaset(dd);
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
%     v=10^-7;
%     u=delta*Z;
%     a=Z;
%     s2=10^-4*ones(nq,1);
    warning off;
    dd
    ss
[~,index8,theta8,theta80,Test,runlength] = thom666(Y1',B00,B,m,sigma,niter,w,sigma9,lambda0,controlimit,sigma0);
% [~,index8,theta8,theta80,Test,runlength] = thom7(Y1',B00,B,m,sigma,controlimit,sigma0,u,a,s2);
TTest(ss,:) = Test;
RRunlength(ss) = runlength;
end
ARL_oc(dd)= mean(RRunlength);
sdrl_oc(dd)=std(RRunlength);
end
delete(gcp('nocreate'))
% save('GenerateVector_vb.mat')
 
% theta = theta8;
% theta0 = theta80;
% index = index8;
% 
% figure('position',[40,100,600,750]);
% %pause(20);
% for i = 2:T
% Yi = Y(:,i);
% Xi1 =B*theta(:,i);
% Xi = Xi1+B00*theta0(:,i);
% subplot(2,2,1);
% imagesc(squeeze(Yi));
% title('Figure 1:Simulation image','Fontsize',8);
% subplot(2,2,2);
% imagesc(Xi);
% title('Figure 2:Fitted image','Fontsize',8);
% subplot(2,2,3);
% imagesc(Xi1);
% title('Figure 3:Inferential anomaly point','Fontsize',8);
% M = zeros(n,1);
% M(index(:,i)) = 1;
% subplot(2,2,4);
% imagesc(M);
% title(int2str(i));
% pause(0.8);
% end


% limit = 0.000134
% for rep = 1:nrep_oc
%     ind =  find(TTest(rep,2:end)>limit,1);
%     if isempty(ind)
%         runlen1(rep) = T;
%     else
%         runlen1(rep) = ind;
%     end
% end
% mean(runlen1)
% std(runlen1)

% for i = 1:n
% num(i) =length(find(index8 == i));
% end
% 
% figure('position',[800,50,20,270]);
% imagesc(theta8);
% 
% figure('position',[800,50,70,500]);
% imagesc(B0)
% 
% plot(Y1(:,80),'LineWidth',4)
% ylim([-0.5,0.6]);
% plot(B*theta8,'r','LineWidth',4)
% ylim([-0.5,0.6]);
% plot(B0*theta80,'g','LineWidth',4)
% ylim([-0.5,0.6]);
% plot(E(:,80),'k','LineWidth',4)
% ylim([-0.5,0.6]);
% plot(Y(:,80),'LineWidth',4)
% ylim([-0.5,0.6]);
% 
% 
% theta8=[zeros(1,9),0.3,zeros(1,7)]';
% theta80=[0.1,0.3,-0.1]';
% plot(Y1(:,80),'LineWidth',4)
% hold on
% plot(B*theta8(:,80),'r','LineWidth',4)
% plot(B0*theta80(:,80),'g','LineWidth',4)
% plot(E(:,80),'k','LineWidth',4)
% plot(Y(:,80),'m','LineWidth',4)
% 
% 
% 
% xlswrite('3.xlsx',Index(:,:,1),'sheet1');
% xlswrite('3.xlsx',Index(:,:,2),'sheet2');
% xlswrite('3.xlsx',Index(:,:,3),'sheet3');
% xlswrite('3.xlsx',Index(:,:,4),'sheet4');
% xlswrite('3.xlsx',Index(:,:,5),'sheet5');
% xlswrite('3.xlsx',Index(:,:,6),'sheet6');
% xlswrite('3.xlsx',Index(:,:,7),'sheet7');
% xlswrite('3.xlsx',Index(:,:,8),'sheet8');
% xlswrite('3.xlsx',Index(:,:,9),'sheet9');
% xlswrite('3.xlsx',Index(:,:,10),'sheet10');