%%% CIV %%%

z=[ones(1000,1),w,s];
x=[ones(1000,1),p,w];

[n,k]=size(x);

sum_zx=zeros(k,k);
sum_zy=zeros(k,1);

for i=1:n
    sum_zx=sum_zx+z(i,:)'*x(i,:);
    sum_zy=sum_zy+z(i,:)'*y(i,1);
end

theta_civ=inv(sum_zx)*sum_zy;

uhat_civ=y-x*theta_civ;

sum_uzz=zeros(k,k);
sum_xz=zeros(k,k);

for i=1:n
    sum_uzz=sum_uzz+uhat_civ(i,1)^2*z(i,:)'*z(i,:)/n;
    sum_xz=sum_xz+x(i,:)'*z(i,:);
end

av_civ=inv(sum_zx)*sum_uzz*inv(sum_xz);
av_civ_m=inv(z'*x)*z'*uhat_civ*uhat_civ'*z*inv(x'*z);
ste_civ=sqrt(diag(av_civ)./n);
sum_zx_m=inv(z'*x);
