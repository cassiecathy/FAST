function basis = fourierbasis(nq,ntime)

point = linspace(-pi,pi,ntime);
norder = ceil(nq/2);
basis = zeros(ntime,nq);
for i = 1:norder
basis(:,2*i-1) = cos(i*point);
if 2*i <= nq
basis(:,2*i) = sin(i*point);
end
end
end