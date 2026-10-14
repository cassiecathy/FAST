clc;clear;

dim = 10;
kmax = 1000;
m0 = 100;
tau = 0;
alphaerr = 0.005;

disttype = 1;
distparams = cell(3);
distparams{1} = zeros(1,dim);

% set cov matrix
   distparams{1} = zeros(1,dim);
    node = 0:dim:dim;
    node = sort(node);
    nodep = diff(node);
    vari = linspace(1,6,5);
    A = cell(1,length(nodep));
    covMat = eye(dim); 
for i = 1:dim
for j = i+1:dim
    covMat(i,j) = 0.5;
    covMat(j,i) = covMat(i,j);
end
end
%
distparams{2} = covMat;
distparams{3} = 5;
distparams1= distparams;

ocdim = 4;
censdim = 5;

distparams1 = distparams;
% distparams1{1}(1:ocdim) = distparams1{1}(1:ocdim) + meanshift;
shift = 0;
g = repmat(1/(dim+1),dim+1,1);
delta = 0.07;
k = 0.05;
lambda = (dim*mean(g(1:dim))-(dim-censdim)*delta)/(1-g(dim+1)^(censdim/dim));
lambda0 = g(dim+1)^(1-censdim/dim);

repNum = 1000;
