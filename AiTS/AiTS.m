function [Tval,runlen,index] = AiTS(data,g,m,m0,k,h,T0)
[k0,p]=size(data);
% ICdata = data(1:m0,:);

OCdata = data(m0+1:end,:);

index = zeros(m,k0-m0);
S1 = zeros(p,1);
S2 = zeros(p,1);
Tval = zeros(k0-m0,1);
omega_sum=zeros(p,1);
runlen=k0;
for i = 1:k0-m0
    omega=g;
    if i==1
        ind = randsample(p,p);
    elseif S1==0
        ind=index(:,i-1);  
    else
        N=S1/sum(S1);
        y=zeros(p,1);
        for j=1:p
            y(j)=gamrnd(N(j),1);
        end
        Z=y/sum(y);
        [val,ind] = sort(Z,'descend');
    end
   curr_ind = ind(1:m);
   index(:,i) = curr_ind;    
   omega(curr_ind)=0;
   currdata = OCdata(i,index(:,i));
   curr_min = find(currdata==max(currdata));
       
   omega(curr_ind(curr_min))=sum(g(curr_ind))*g(curr_ind(curr_min))/sum(g(curr_ind(curr_min)));

   omega_sum=omega_sum+omega; 
   theta=(g+omega_sum)/(i+1);
   
   eta=theta;
   eta(curr_ind)=0;
   eta(curr_ind(curr_min))=sum(theta(curr_ind));
   
   Ct = (S1-S2+eta-g)'*diag(1./(S2+g))*(S1-S2+eta-g);
   if Ct < k
       S1 = zeros(p,1);
       S2 = zeros(p,1);
   else
       S1 = (S1+eta).*(Ct-k)./Ct;
       S2 = (S2+g).*(Ct-k)./Ct;
   end
   Tval(i) = (S1-S2)'*diag(1./S2)*(S1-S2);
   if Tval(i)>h && i>T0
       runlen = i-T0;
       break;
   end 
    
end

end
 