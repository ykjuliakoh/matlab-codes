sigma_sq=exp(2*x(:,2)+4*x(:,3));

[n,k]=size(x);

sum_xx=zeros(k,k);
sum_xy=zeros(k,1);

for i=1:n
    sum_xx=sum_xx+x(i,:)'*x(i,:)/(sigma_sq(i,1)*n);
    sum_xy=sum_xy+x(i,:)'*y(i,1)/(sigma_sq(i,1)*n);
end

theta_gls=inv(sum_xx)*sum_xy;

uhat=y-x*theta_gls;

sum_uxx=zeros(k,k);

for i=1:n
    sum_uxx=sum_uxx+uhat(i,1)^2*x(i,:)'*x(i,:)/(sigma_sq(i,1)*n);
end

av_gls=inv(sum_xx)*sum_uxx*inv(sum_xx);
ste_gls=sqrt(diag(av_gls)./n);
