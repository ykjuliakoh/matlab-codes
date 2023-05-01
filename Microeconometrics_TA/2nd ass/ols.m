%%% OLS %%%

[n,k]=size(x);

sum_xx=zeros(k,k);
sum_xy=zeros(k,1);

for i=1:n
    sum_xx=sum_xx+x(i,:)'*x(i,:)/n;
    sum_xy=sum_xy+x(i,:)'*y(i,1)/n;
end

theta_ols=inv(sum_xx)*sum_xy;

uhat=y-x*theta_ols;

sum_uxx=zeros(k,k);

for i=1:n
    sum_uxx=sum_uxx+uhat(i,1)^2*x(i,:)'*x(i,:)/n;
end

sum_uxx_v=sum_uxx;
Sigma=diag(diag(uhat*uhat'));
sum_uxx_m=x'*Sigma*x;
av_ols=inv(sum_xx)*sum_uxx*inv(sum_xx);
av_ols_m=inv(x'*x)*sum_uxx_m*inv(x'*x);
ste_ols=sqrt(diag(av_ols)./n);
