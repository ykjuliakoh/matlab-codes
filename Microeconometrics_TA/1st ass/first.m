x=load('x.dat');
y=load('y.dat');

[N,K]=size(x); %Size of x, N=observations and 
%K=numbers of explanatory variables

% (a) Determine OLS estimator beta_hat.
sum_xx=zeros(K,K);
sum_xy=zeros(K,1);

for i=1:N
    sum_xx=sum_xx+x(i,:)'*x(i,:);
    sum_xy=sum_xy+x(i,:)'*y(i,1);
end

sum_xx=sum_xx/N;
sum_xy=sum_xy/N;

beta_hat=inv(sum_xx)*sum_xy;

% (b), (c) Determine standard erros of beta_hat.
u_hat=y-x*beta_hat;

% Homoskedasticity
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

AV_hetero=inv(A_hetero)*B_hetero*inv(A_hetero);
SE_hetero=sqrt(diag(AV_hetero)/N);

% (c) partial coefficient
x_temp=[];
x_temp(:,1:2)=x(:,1:2);
x_temp(:,3)=x(:,4);
x2=x(:,3);

beta_temp=inv(x_temp'*x_temp)*(x_temp'*y);
resid_y=y-x_temp*beta_temp;

beta_temp2=inv(x_temp'*x_temp)*(x_temp'*x2);
resid_x2=x2-x_temp*beta_temp2;

beta_tilde=inv(resid_x2'*resid_x2)*(resid_x2'*resid_y);
chk=beta_tilde-beta_hat(3,1);

% (e) t-test
t_test_homo=beta_hat./SE_homo;
t_test_hetero=beta_hat./SE_hetero;

% (f) general t-test eg1
q_f=[0;1;1;1];

vq_f_homo=q_f'*AV_homo*q_f;
t_q_f_homo=sqrt(N)*(beta_hat(2,1)+beta_hat(3,1)+beta_hat(4,1))/vq_f_homo;

vq_f_hetero=q_f'*AV_hetero*q_f;
t_q_f_hetero=sqrt(N)*(beta_hat(2,1)*beta_hat(3,1)*beta_hat(4,1))/vq_f_hetero;

% (g) general t-test eg2
q_g=[0;beta_hat(4,1);-2*beta_hat(3,1);beta_hat(2,1)];

vq_g_homo=q_g'*AV_homo*q_g;
t_q_g_homo=sqrt(N)*(beta_hat(2,1)*beta_hat(4,1)-beta_hat(3,1)^2+1)/vq_g_homo;

vq_g_hetero=q_g'*AV_hetero*q_g;
t_q_g_hetero=sqrt(N)*(beta_hat(2,1)*beta_hat(4,1)-beta_hat(3,1)^2+1)/vq_g_hetero;

% (h) F-test: J=3 (# of restriction=3)
q_h=[beta_hat(2,1);beta_hat(3,1);beta_hat(4,1)];
Q_h=[0,1,0,0;0,0,1,0;0,0,0,1];

F_q_h_homo=N*q_h'*inv(Q_h*AV_homo*Q_h')*q_h;
F_q_h_hetero=N*q_h'*inv(Q_h*AV_hetero*Q_h')*q_h;

% (i) General F-test: J=2
q_i=[beta(2,1)-0.5;beta_hat(3,1)+beta_hat(4,1)+0.5];
Q_i=[0,1,0,0;0,0,1,1];

F_q_i_homo=N*q_i'*inv(Q_i*AV_homo*Q_i')*q_i;
F_q_i_hetero=N*q_i'*inv(Q_i*AV_hetero*Q_i')*q_i;

% (j) General F-test: J=2
q_j=[beta_hat(3,1)-2;beta_hat(2,1)*beta_hat(4,1)-1];
Q_j=[0,0,1,0;0,beta_hat(4,1),0,beta_hat(2,1)];

F_q_j_homo=N*q_j'*inv(Q_j*AV_homo*Q_j')*q_j;
F_q_j_hetero=N*q_j'*inv(Q_j*AV_hetero*Q_j')*q_j;
   
% (k) R^2
SSR=zeros(1,1);
SST=zeros(1,1);

for i=1:N
    SSR=SSR+u_hat(i,1)^2;
    SST=SST+(y(i,1)-mean(y))^2;
end

SSR=SSR/N;
SST=SST/N;

R_sq=1-SSR/SST;

