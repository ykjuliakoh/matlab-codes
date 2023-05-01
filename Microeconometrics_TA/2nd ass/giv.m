x=[ones(1000,1),p,w];
z=[ones(1000,1),w,s,t];

[n,k]=size(x);
[n,dz]=size(z);

sum_xz=zeros(k,dz);
sum_zz=zeros(dz,dz);
sum_zx=zeros(dz,k);
sum_zy=zeros(dz,1);

for i=1:n
    sum_xz=sum_xz+x(i,:)'*z(i,:)/n;
    sum_zz=sum_zz+z(i,:)'*z(i,:)/n;
    sum_zx=sum_zx+z(i,:)'*x(i,:)/n;
    sum_zy=sum_zy+z(i,:)'*y(i,1)/n;
end

theta_giv=inv(sum_xz*inv(sum_zz)*sum_zx)*sum_xz*inv(sum_zz)*sum_zy;

uhat_giv=y-x*theta_giv;

sum_uzz=zeros(dz,dz);

for i=1:n
    sum_uzz=sum_uzz+uhat_giv(i,1)^2*z(i,:)'*z(i,:)/n;
end

av_giv=inv(sum_xz*inv(sum_zz)*sum_zx)*sum_xz*inv(sum_zz)*sum_uzz*inv(sum_zz)*sum_zx*inv(sum_xz*inv(sum_zz)*sum_zx);
ste_giv=sqrt(diag(av_giv)./n);
