run('EsettingstdNA.m');
tic;

ocdim = [1,2,3,4];
Tval=zeros(repNum, kmax-m0);

% % h = 4.7; % t3
% h = 9.76; % norm delta = 0.01
% h = 5.94; % norm, p = 12, covmat = eye(dim)
% h = 6.37;  % norm, p = 12, delta = 0.08, covMat = 0.5, 
% h = 5.5; %norm
deltaset = 0.25*[0,1,2,3,4,5,6,7,8]'*[1,-1,1,-1];
deltaset= abs(deltaset);
mumin = 0.25;
h = inf;
for dd = 1:size(deltaset,1)
    shift = deltaset(dd,:);
    distparams1{1} = distparams{1};
    distparams1{1}(ocdim)=distparams{1}(ocdim)+shift;
    parfor i = 1 : repNum
        tic;
        data = zeros(kmax, dim);
        data(1:m0+tau,:)=RandomSample((m0+tau), disttype, distparams);
        data((m0+tau+1):kmax,:) = RandomSample((kmax-tau-m0), disttype, distparams1);
       [Tval(i,:),runlen(dd,i)] = naCUSUM(data,censdim,m0,delta,k,h,lambda,lambda0,g);
    end
    arl(dd) = mean(runlen(dd,:));
    stdrl(dd) = std(runlen(dd,:));
end

mean(runlen,2)

save(sprintf('%s_Total%f.mat', filename, mod(now,1)));
toc;
fclose('all');

h =29.9
for i = 1:repNum
    ind = find(Tval(i,:)>h,1);
    if ind< kmax-m0
        runlen(1,i) = ind;
    else
        runlen(1,i) = kmax-m0;
    end
end
mean(runlen(1,:))
std(runlen(1,:))
