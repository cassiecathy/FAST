function [Tval,runlen,index] = naCUSUM(data,censdim,m0,delta,k,h,lambda,lambda0,g)

[k0,p] = size(data);
ICdata = data(1:m0,:);
ICdata1 = [ICdata,zeros(m0,1)];

OCdata = data(m0+1:end,:);
OCdata = [OCdata,zeros(k0-m0,1)];
S1 = zeros(p+1,1);
SS1 = zeros(p+1,k0-m0);
SS2 = zeros(p+1,k0-m0);
xxi = zeros(p+1,k0-m0);
S2 = zeros(p+1,1);
curr_indd = zeros(censdim,k0-m0);
Tval = zeros(1,k0-m0);
for i = 1:k0-m0
   xi = delta.*ones(p+1,1);
   if i ==1
    ind = randsample(p+1,p+1);
   else
   [val, ind] = sort(S1,'descend');
   end
   curr_ind = ind(1:censdim);
   curr_indd(:,i) = curr_ind;
   xi(curr_ind) = 0;
   index(i,:) = curr_ind;
   currdata = [OCdata(i,index(i,:)),0];
   curr_min = find(currdata==min(currdata));
   if curr_min< censdim+1
   xi(curr_ind(curr_min)) = lambda;
   xi(p+1) = 0;
   else
   xi(p+1) = lambda0;    
   end
   xxi(:,i) = xi;
   Ct = (S1-S2+xi-g)'*diag(1./(S2+g))*(S1-S2+xi-g);
   if Ct < k
       S1 = zeros(p+1,1);
       S2 = zeros(p+1,1);
   else
       S1 = (S1+xi).*(Ct-k)./Ct;
       S2 = (S2+g).*(Ct-k)./Ct;
   end
   Tval(i) = (S1-S2)'*diag(1./S2)*(S1-S2);
   if Tval(i) >h && i>1
       break;
   end
   SS1(:,i) = S1;
   SS2(:,i) = S2;
end
runlen = i;
end
 