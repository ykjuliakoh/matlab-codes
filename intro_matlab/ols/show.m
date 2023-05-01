clc
clear

x=load('x.dat');
y=load('y.dat');

[N,K]=size(x);

% Determine OLS estimator beta_hat.
sum_xx=zeros(K,K);
sum_xy=zeros(K,1);

for i=1:N
    sum_xx=sum_xx+x(i,:)'*x(i,:);
    sum_xy=sum_xy+x(i,:)'*y(i,1);
end

beta_hat=inv(sum_xx)*sum_xy;%individual by individual calculation

b_hat=inv(x'*x)*x'*y; %easier version

b_hat2 = (x'*x)\(x'*y);


%Determine standard erros of beta_hat.
u_hat=y-x*b_hat;

% Homoskedasticity
% individual by individual summation

sum_homo=zeros(K,K);
sum_sigmasq=zeros(1,1);

for i=1:N
    sum_homo=sum_homo+x(i,:)'*x(i,:);
    sum_sigmasq=sum_sigmasq+u_hat(i,1)^2;
end

sum_homo=sum_homo/N;
sum_sigmasq=sum_sigmasq/N;

AV_homo=sum_sigmasq*inv(sum_homo);
SE_homo=sqrt(diag(AV_homo)/N);

% Heteroskedasticity
A_hetero=zeros(K,K);
B_hetero=zeros(K,K);

for i=1:N
    A_hetero=A_hetero+x(i,:)'*x(i,:);
    B_hetero=B_hetero+u_hat(i,1)^2*x(i,:)'*x(i,:);
end

A_hetero=A_hetero/N;
B_hetero=B_hetero/N;

A_het = (x'*x)/N;
B_hat = ((u_hat.*x)'*(u_hat.*x))/N;

AV_hetero=inv(A_hetero)*B_hetero*inv(A_hetero);
SE_hetero=sqrt(diag(AV_hetero)/N);

% R^2
SSR=zeros(1,1);
SST=zeros(1,1);

for i=1:N
    SSR=SSR+u_hat(i,1)^2;
    SST=SST+(y(i,1)-mean(y))^2;
end

SSR2 = u_hat'*u_hat;
SSR=SSR/N;
SST=SST/N;

R_sq=1-SSR/SST;
