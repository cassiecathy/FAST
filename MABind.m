function [runlen,Test,count1,Q,meanU1,index] = MABind(Yic,Yoc,m,limit,tau0)
m0 = 0;
lambda = 0.9;

data = Yoc';
[T,dim] = size(data);
normal_mean = mean(Yic,2);
covMat = 1/(T-1) * (Yic - repmat(normal_mean,1,T)) * (Yic - repmat(normal_mean,1,T))';
loss = @(t) sqrt(max(0,log((1-(1-lambda).^(t))./lambda)));

meanU0 = zeros(1,dim);
invcovMat = inv(covMat + 1e-4*eye(dim)); %#ok<NASGU>

Q = zeros(T-m0,dim);
count = zeros(T-m0,dim);
runlen = T-m0-tau0;

index = zeros(m,T-m0);
meanU = zeros(T-m0,dim);
sumU = zeros(T-m0,dim);
Omega = cell(T-m0,1);
ww = 1-lambda;
Test = zeros(1,T-m0);

for i = 1:T-m0
    if i == 1
        % First step: also use partial observation only
        index(:,1) = randsample(dim,m,false)';
        tmp = sort(index(:,1));

        obsM = data(m0+1,tmp);

        dataM = zeros(1,dim);
        dataM(tmp) = obsM;

        Z = eye(dim);
        Z = Z(tmp,:);

        Omega{1} = Z' * inv(covMat(tmp,tmp)) * Z;
        sumU(1,:) = meanU0 + dataM(tmp) * inv(covMat(tmp,tmp)) * Z;
        meanU(1,:) = sumU(1,:) * inv(Omega{1} + 1e-4*eye(dim));

        count(1,tmp) = 1;
        ww(1) = 1;
    end

    if i > 1
        ww(1:i-1) = ww(1:i-1) .* (1-lambda);
        ww(i) = 1;

        [~,ccc] = sort(Q(i-1,:),'descend');
        index(:,i) = ccc(1:m)';

        count(i,:) = count(i-1,:);
        count(i,index(:,i)) = count(i,index(:,i)) + 1;

        obsM = data(i+m0,index(:,i));
        non_ind = 1:dim;
        non_ind(index(:,i)) = [];

        dataM = zeros(1,dim);
        dataM(index(:,i)) = obsM;

        Z = eye(dim);
        Z(non_ind,:) = [];
        tmp = sort(index(:,i));

        Omega{i} = Omega{i-1} * (1-lambda) + Z' * inv(covMat(tmp,tmp)) * Z;
        sumU(i,:) = sumU(i-1,:) * (1-lambda) + dataM(tmp) * inv(covMat(tmp,tmp)) * Z;
        meanU(i,:) = sumU(i,:) * inv(Omega{i} + 1e-4*eye(dim));
    end

    Test(i) = 1 ./ sum(ww(1:i).^2) * sum(ww(1:i)) * meanU(i,:) * Omega{i} * meanU(i,:)';


    if Test(i) > limit && i > tau0
        runlen = i - tau0;
        break;
    end

    Q(i,:) = (meanU(i,:) + sign(meanU(i,:)) .* loss(i) .* sqrt(diag(inv(Omega{i}+1e-4*eye(dim))))').^2;
end

count1 = count(runlen+tau0,:);
meanU1 = meanU(tau0+runlen,:);

end