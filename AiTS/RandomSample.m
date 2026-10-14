function sample = RandomSample(n, type, params)
% generating random number specified by the type and parameters
% n is the number of samples required
% type is the distribution type
% params is the cell array containing corresponding paramemters
% output: sample is the matrix n-by-dim

% multivariate normal distribution 
if type == 1
    sample = mvnrnd(params{1}, params{2}, n);
    return;
end

% multivariate t-distribution 
if type == 2
    dim = length(params{1});
    sample = mvnrnd(zeros(n,dim), params{2});
    denom = chi2rnd(params{3}, n, 1); % params{3} is the degree of freedom
    denom = sqrt(denom/params{3});
    sample = sample./repmat(denom, 1, dim);
    sample = sample+repmat(params{1}, n, 1);
    return;
end

% multivariate gamma distribution 
if type == 3
    dim=length(params{1})
    sample = zeros(n, dim);
    for i = 1 : n
        tmp = mvnrnd(zeros(1,dim), params{2}, params{3});
        sample(i,:) = sum(tmp.^2)/2+params{1};
    end
    return;
end
end