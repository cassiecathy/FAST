clc;
clear;

n = 100; % Image size n*n
T = 1000; % The Total time length
bd = 4; % How smooth the defect is (Order of Spline)
% k0 = 20; % scale of the background
% bd0 = 3;
k1  = 4; % scale of the defect
ndefect = 1;  % How many defect in the system
delta =0; % Defect magnitude
sigma = 0.05; % The noise level
Tin =0; % When the change occur
%Y = zeros(n,n,T); % Generated Videos
%cutoff = 0.5;
nq0 =3;


repNum =1000;
tau=0;
m =30; % budget of pixels that can be observed for each image
lambda=0.1;
limit=1500;
m0=0;

% for i = 1:p
% for j = i+1:p
%     covMat(i,j) = 0.5^(abs(i-j));
%     covMat(j,i) = covMat(i,j);
% end
% end

loss = @(t) sqrt(max(0,log((1-(1-lambda).^t)./lambda)));
runlen=zeros(1,repNum);
Test = zeros(repNum,T-m0);
count1 = zeros(repNum,n);
meanU = count1;
deltaset=[0.5,0.8,1.0,1.5,2.0];
arl_mab=zeros(1,size(deltaset,2));
sdrl_mab=zeros(1,size(deltaset,2));

parpool('local',8);
for dd=1:size(deltaset,2)
    delta=deltaset(dd);
for rep =  1:repNum
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
    covMat=eye(n)*sigma^2+sigma0.*B00*eye(size(B00,2))*B00';
    
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
    Y=(Y-repmat(normal_mean,1,T))./normal_std;
    dd
    rep
   %tic;
   [count1,Q,meanU,runlen,Test,index] = MABind(Y',m0,m,covMat,loss,lambda,limit,tau);
   %toc;
TTest(rep,:) = Test;
rrunlen(rep,:) = runlen;
% Index(:,:,rep)=index;
end
arl_mab(dd)= mean(rrunlen);
sdrl_mab(dd)= std(rrunlen);
end
delete(gcp('nocreate'))
% save('MABmain.mat')

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

% figure('position',[800,500,1000,400]);
% for i=2:T
%     subplot(1,2,1);
%     imagesc(squeeze(Y(:,:,i)));
%     M=zeros(n,n);
%     M(index(:,i))=1;
%     subplot(1,2,2);
%     imagesc(M);
%     title(int2str(i));
%     pause(0.5);
% end

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