function      [index_ini,index,theta,theta0,Test,runlength] = thom666(data,trueA0,trueA,m,sigma,niter,w,sigma9,lambda0,controlimit,sigma0,T0)
[n,p] = size(data);
[~,q0] = size(trueA0);
[~,q] = size(trueA);

index = zeros(m,n);
B0=zeros(m,q0,n);
B1=zeros(m,q,n);
X0=zeros(m,n);

XB=zeros(q,n);
BB=zeros(q,n);
BBjk=zeros(q,q,n);
theta=zeros(q,n);
theta0=zeros(q0,n);

Test = zeros(1,n);
runlength = n;
time=zeros(n,1);

u_out=zeros(n,q);
a_out=zeros(n,q);
s2_out=zeros(n,q);


%% for the first n=1
index_ini =1:p;
B0_ini = trueA0(index_ini,:);
B1_ini = trueA(index_ini,:);     
BB(:,1)=diag(B1_ini'*B1_ini);
BBjk(:,:,1)=B1_ini'*B1_ini;   
s2=1./(BB(:,1)/sigma^2+1./sigma9.^2);
a=0.5*ones(q,1);
u=zeros(q,1);
 
X0_ini = data(1,index_ini)';
mu_tmp = zeros(q,1);
inv_B0B0 = inv(B0_ini'*B0_ini+1e-6.*eye(q0));
theta0(:,1) = inv_B0B0*(B0_ini'*(X0_ini-B1_ini*mu_tmp));
X1_ini = X0_ini - B0_ini*theta0(:,1);
for h=1:niter    
    XB(:,1) = (X1_ini'*B1_ini)';   
     ss=zeros(q,1);
    for j=1:q
    u(j)=(XB(j,1)-(BBjk(j,:,1)*(a.*u)-BBjk(j,j,1)*a(j)*u(j)))*s2(j)/sigma^2;

    ss(j)=log(w(j)/(1-w(j)))+1/2-(s2(j)+u(j)^2)/(2*sigma9(j)^2)+log(sqrt(s2(j))/sigma9(j))+u(j)^2/s2(j)-BB(j,1)*(u(j)^2+s2(j))/(2*sigma^2);
%     ss(j)=log(w(j)/(1-w(j)))-u(j)^2/(2*sigma9(j)^2)+u(j)^2/s2(j)-BB(j,1)*(u(j)^2+s2(j))/(2*sigma^2);
    a(j)=exp(ss(j))/(1+exp(ss(j)));

    if isnan(a(j))==1
        a(j)=1;
    end
    end

     mu_tmp = u.*a;
     sigma_b(:,:,1)=(B0_ini'*(sigma^(-2))*B0_ini+sigma0^(-1).*eye(q0))^(-1);
     mu_theta0(:,1)=sigma_b(:,:,1)*B0_ini'*(sigma^(-2))*(X0_ini-B1_ini*mu_tmp);
     X1_ini = X0_ini - B0_ini*mu_theta0(:,1);
%      theta0(:,1) = inv_B0B0*(B0_ini'*(X0_ini-B1_ini*mu_tmp));
%      X1_ini = X0_ini - B0_ini*theta0(:,1);
end


r=zeros(q,1);
for j=1:q
    r(j)=random('Binomial',1,a(j));
    if r(j)==1
        theta(j,1)=normrnd(u(j),s2(j)^0.5);
    else
        theta(j,1)=0;
%         theta(j,1)=normrnd(0,(v*s2(j))^0.5);
    end
end

A=zeros(q,q);
for j=1:q
    for k=1:q
        A(j,k)=a(j)*a(k);
    end
    A(j,j)=A(j,j)^(1/2);
end
theta_a = u.*a;
u_out(1,:)=u;
a_out(1,:)=a;
s2_out(1,:)=s2;

E=normrnd(0,sigma,[1,p]);
X=trueA*theta(:,1)+E';
thom=zeros(p,1);
for k=1:p
thom(k)=2*X(k)*(trueA(k,:)*theta_a)-trueA(k,:)*(u*u'.*A)*trueA(k,:)';
end
shuffle = randsample(1:p,p,0);
[~,ind] = sort(thom(shuffle),'descend');
 index(:,2)=shuffle(ind(1:m));
Test(1)=2*X1_ini'*squeeze(B1_ini)*theta_a-u'*(squeeze(B1_ini)'*squeeze(B1_ini).*A)*u;

%% for later n = 2, 3,....
for i = 2:n
t0=clock;
     B0(:,:,i)= trueA0(index(:,i),:);
     B1(:,:,i)= trueA(index(:,i),:);
 %    X0(:,i)= data(i,index(:,i))';
     
        BB(:,i)=BB(:,i-1)*lambda0+diag(B1(:,:,i)'*B1(:,:,i));
        BBjk(:,:,i)=BBjk(:,:,i-1)*lambda0+B1(:,:,i)'*B1(:,:,i);

 s2=1./(BB(:,i)/sigma^2+1./sigma9.^2);
 a=0.5*ones(q,1);
 u=zeros(q,1);
 
X0(:,i)= data(i,index(:,i))';
mu_tmp = zeros(q,1);
inv_B0B0 = inv(B0(:,:,i)'*B0(:,:,i)+1e-6.*eye(q0));
theta0(:,i) = inv_B0B0*(B0(:,:,i)'*(X0(:,i)-B1(:,:,i)*mu_tmp));
X1(:,i)= X0(:,i) - B0(:,:,i)*theta0(:,i);
for h=1:niter    
    if i ==1
      XB(:,i) = (X1(:,i)'*B1(:,:,i))';   
    else
      XB(:,i)=XB(:,i-1)*lambda0+(X1(:,i)'*B1(:,:,i))';          
    end
     ss=zeros(q,1);

    for j=1:q
    u(j)=(XB(j,i)-(BBjk(j,:,i)*(a.*u)-BBjk(j,j,i)*a(j)*u(j)))*s2(j)/sigma^2;
% fid=fopen('test3.txt','at');
% fprintf(fid,'%f ',u(j));
% fclose(fid);

     ss(j)=log(w(j)/(1-w(j)))+1/2-(s2(j)+u(j)^2)/(2*sigma9(j)^2)+log(sqrt(s2(j))/sigma9(j))+u(j)^2/s2(j)-BB(j,i)*(u(j)^2+s2(j))/(2*sigma^2);
%     ss(j)=log(w(j)/(1-w(j)))-u(j)^2/(2*sigma9(j)^2)+u(j)^2/s2(j)-BB(j,i)*(u(j)^2+s2(j))/(2*sigma^2);

    a(j)=exp(ss(j))/(1+exp(ss(j)));

    if isnan(a(j))==1
        a(j)=1;
    end
    end

     mu_tmp = u.*a;
     sigma_b(:,:,i)=(B0(:,:,i)'*(sigma^(-2).*eye(m))*B0(:,:,i)+sigma0^(-1).*eye(q0))^(-1);
     mu_theta0(:,i)=sigma_b(:,:,i)*B0(:,:,i)'*(sigma^(-2).*eye(m))*(X0(:,i)-B1(:,:,i)*mu_tmp);
     X1(:,i) = X0(:,i) - B0(:,:,i)*mu_theta0(:,i);
%      theta0(:,i) = inv_B0B0*(B0(:,:,i)'*(X0(:,i)-B1(:,:,i)*mu_tmp));
%      X1(:,i) = X0(:,i) - B0(:,:,i)*theta0(:,i);
%      fid=fopen('test3.txt','at');
%      fprintf(fid,'\n');
%      fprintf(fid,'\n');
%      fclose(fid);
end


r=zeros(q,1);
for j=1:q
    r(j)=random('Binomial',1,a(j));
    if r(j)==1
        theta(j,i)=normrnd(u(j),s2(j)^0.5);
    else
      theta(j,i)=0;
%         theta(j,i)=normrnd(0,(v*s2(j))^0.5);
    end
end
for j=1:q
    for k=1:q
        A(j,k)=a(j)*a(k);
    end
    A(j,j)=A(j,j)^(1/2);
end

theta_a = u.*a;
u_out(i,:)=u;
a_out(i,:)=a;
s2_out(i,:)=s2;


E=normrnd(0,sigma,[1,p]);
X=trueA*theta(:,i)+E';

thom=zeros(p,1);
for k=1:p
thom(k)=2*X(k)*(trueA(k,:)*theta_a)-trueA(k,:)*(u*u'.*A)*trueA(k,:)';
end
shuffle = randsample(1:p,p,0);
[~,ind] = sort(thom(shuffle),'descend');
index(:,i+1)=shuffle(ind(1:m));
% index(:,i+1)=randsample(p,m);
 
Test(i)=2*X1(:,i)'*squeeze(B1(:,:,i))*theta_a-u'*(squeeze(B1(:,:,i))'*squeeze(B1(:,:,i)).*A)*u;
if Test(i)>controlimit && i>T0
runlength =i-T0;
break;
end
t1=clock;
time(i)=etime(t1,t0);
end
end






