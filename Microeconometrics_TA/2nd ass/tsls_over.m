%%% TSLS of overidentification %%%

%%% 1st stage %%%
x=[ones(1000,1),p,w];
z=[ones(1000,1),w,s,t];

[n,k]=size(x);
[n,dz]=size(z);

sum_zz=zeros(dz,dz);
sum_zp=zeros(dz,1);

for i=1:n
    sum_zz=sum_zz+z(i,:)'*z(i,:)/n;
    sum_zp=sum_zp+z(i,:)'*p(i,1)/n;
end

p_ols=inv(sum_zz)*sum_zp;

phat=z*p_ols;

%%% 2nd stage %%%

x=[ones(1000,1),phat,w];

sum_xx=zeros(k,k);
sum_xy=zeros(k,1);

for i=1:n
    sum_xx=sum_xx+x(i,:)'*x(i,:)/n;
    sum_xy=sum_xy+x(i,:)'*y(i,1)/n;
end

theta_2sls=inv(sum_xx)*sum_xy;

uhat=y-x*theta_2sls;

sum_uxx=zeros(k,k);

for i=1:n
    sum_uxx=sum_uxx+uhat(i,1)^2*x(i,:)'*x(i,:)/n;
end

av_2sls=inv(sum_xx)*sum_uxx*inv(sum_xx);
ste_2sls=sqrt(diag(av_2sls)./n);


