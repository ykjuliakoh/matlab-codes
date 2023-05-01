z=[ones(1000,1),w,s];

[n,k]=size(z);

%%% 1st stage %%%
sum_zz=zeros(k,k);
sum_zp=zeros(k,1);

for i=1:n
    sum_zz=sum_zz+z(i,:)'*z(i,:)/n;
    sum_zp=sum_zp+z(i,:)'*p(i,1)/n;
end

p_ols=inv(sum_zz)*sum_zp;

phat=z*p_ols;

%%% 2nd stage %%%
x=[ones(1000,1),phat,w];

[n,k]=size(x);

sum_xx=zeros(k,k);
sum_xy=zeros(k,1);

for i=1:n
    sum_xx=sum_xx+x(i,:)'*x(i,:)/n;
    sum_xy=sum_xy+x(i,:)'*y(i,1)/n;
end

theta_2sls=inv(sum_xx)*sum_xy;

uhat_2sls=y-x*theta_2sls;

sum_uxx=zeros(k,k);

for i=1:n
    sum_uxx=sum_uxx+uhat_2sls(i,1)^2*x(i,:)'*x(i,:)/n;
end

av_2sls=inv(sum_xx)*sum_uxx*inv(sum_xx);
ste_2sls=sqrt(diag(av_2sls)./n);
