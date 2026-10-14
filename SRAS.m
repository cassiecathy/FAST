function [runlen,Q,index] = SRAS(X_nor,X,cl,m,tau0)
h=1.1;
lambda = 0.1;
delta = 0.1;
[p,T]=size(X);
[~,m0]=size(X_nor);
index = zeros(m,T);
uindex = zeros(p-m,T);
X_tilde=zeros(p,1);

Sigma0=cov(X_nor');
[V,D]=eig(Sigma0);
M0=V*diag(diag(D).^(-1/2));
R0_nor=zeros(m0,1);
for j=1:m0
    for i=1:m0
        if norm(M0*(X_nor(:,j)-X_nor(:,i)))==0
            R0_nor(j)=0;
        else
            R0_nor(j)=1;
        end
    end
end
xi=mean(R0_nor);
R=zeros(p,T);
v=zeros(p,T);
Q=zeros(T,1);
L=zeros(p,T);
z=zeros(p,T);
runlen=T-tau0;


for t = 1:T
    if t==1
        ind = randsample(p,p); 
    else
        [val,ind] = sort(z(:,t-1),'descend');
    end
    index(:,t) =ind(1:m);
    uindex(:,t)=setdiff(1:p,index(:,t));
    X_o=X(index(:,t),t);
    X_u=X(uindex(:,t),t);
    X_tilde(index(:,t))=X_o;
    KX=zeros(p-m,m0);
    K=zeros(m0,1);
    for i=1:m0
       KX(:,i)=exp(-sum((X_o-X_nor(index(:,t),i)).^2)/h^2/2)*X_nor(uindex(:,t),i);
       K(i)=exp(-sum((X_o-X_nor(index(:,t),i)).^2)/h^2/2);
    end
    X_tilde(uindex(:,t))=sum(KX,2)/sum(K);
    x_sta=zeros(p,m0);
    R0=zeros(p,m0);
    for i=1:m0
        x_sta(:,i)=M0*(X_tilde-X_nor(:,i));
        if norm(x_sta(:,i))==0
            R0(:,i)=zeros(p,1);
        else
            R0(:,i)=x_sta(:,i)/norm(x_sta(:,i));
        end
    end
    
    R(:,t)=mean(R0,2);
    if t==1
        v(:,t)=lambda*R(:,t);
    else
        v(:,t)=(1-lambda)*v(:,t-1)+lambda*R(:,t);
    end
    Q(t)=(2-lambda)*p/(lambda*xi)*norm(v(:,t))^2;
    if Q(t)>cl && t>tau0
       runlen = t-tau0;
       break;
    end 
   
    L(index(:,t),t)=R(index(:,t),t).^2;
    L(uindex(:,t),t)=R(uindex(:,t),t).^2+delta;
    if t==1
        z(:,t)=lambda*L(:,t);    
    else
        z(:,t)=(1-lambda)*z(:,t-1)+lambda*L(:,t);
    end
end

end