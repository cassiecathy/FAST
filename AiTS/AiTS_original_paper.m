clc;
clear;

p=50;
T_nor=1000;
Y_nor=zeros(p,T_nor);
for i=1:p/5
    Y_nor(i,:)=normrnd(0,1,1,T_nor);
end

for i=p/5+1:p/5*2
    Y_nor(i,:)=trnd(3,1,T_nor);
end

for i=p/5*2+1:p/5*3
    Y_nor(i,:)=exprnd(1,1,T_nor);
end

for i=p/5*3+1:p/5*4
    Y_nor(i,:)=lognrnd(1,0.5,1,T_nor);
end

for i=p/5*4+1:p
    Y_nor(i,:)=poissrnd(10,1,T_nor);
end
normal_mean=mean(Y_nor,2);
normal_std=std(Y_nor,0,2);
Y_nor1=(Y_nor-normal_mean)./normal_std;
[M,I]=min(Y_nor1);
B = tabulate(I);
g=B(:,3);


T=1000;
Y=zeros(p,T);
delta=2;
for i=1:p/5
    Y(i,:)=normrnd(0-delta,1,1,T);
end

for i=p/5+1:p/5*2
    Y(i,:)=trnd(3,1,T);
end

for i=p/5*2+1:p/5*3
    Y(i,:)=exprnd(1,1,T);
end

for i=p/5*3+1:p/5*4
    Y(i,:)=lognrnd(1,0.5,1,T);
end

for i=p/5*4+1:p
    Y(i,:)=poissrnd(10,1,T);
end
Y1=(Y-normal_mean)./normal_std;

[test,runlen,index] = AiTS(Y',g,m,m0,k,h,T0);
