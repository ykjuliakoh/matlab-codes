global x y x1 x2 N K

randn ('seed',4);
xu=randn(1000,3);
x1=xu(:,1); x2=xu(:,2);
u=xu(:,3);
%beta0=[1;-2;3];
beta0=[0.1;-0.2;0.3];
x=[ones(1000,1),x1,x2];
[N,K]=size(x);
y=[];
for j=1:1000
    y(j,1)=exp(x(j,:)*beta0) + u(j,1);
end

% OLS %

b_ols=inv(x'*x)*x'*y;   % absurd value

% NLLS %
% using function "grad_nlls.m" and "hess_nlls.m"
theta_ini=[1,1,1]';
% Or try anther initial value
% theta_ini=b_ols;

%test_g=grad_nlls(theta);
%test_h=hess_nlls(theta);

theta_old=theta_ini;
theta_new=theta_old-inv(hess_nlls(theta_old)) * grad_nlls(theta_old);

iter=1;
dist= norm(theta_new-theta_old);

while dist>0.00001
  theta_old=theta_new -inv(hess_nlls(theta_new))*grad_nlls(theta_new);
  iter=iter+1;
  dist= norm(theta_new-theta_old);
  theta_new=theta_old;
   if iter>2000
       break
   end

end
theta_nl=theta_new;


%%%%%%%%%%%%%%% Simplex Method %%%%%%%%%%%%
theta0=[1,1,1];
format long; warning off all; optimset('Tolfun',1e-20,'Tolx',1e-18,'MaxFunEvals',300*length(theta0));
[theta, obj, flag]=fminsearch(@sample_obj_nlls,theta0); 
result_simplex_method=[theta, obj,flag];
theta_nl_simplex=theta;


%%%%%%%%%%%%%%%%% Inference %%%%%%%%%%%%%%%%
u_hat=[];
A_hat=zeros(K,K);
B_hat=zeros(K,K);
for j=1:1000
    x_j=x(j,:)';
    u_hat_j=y(j,1)-exp(x_j'*theta_nl);
    A_hat=A_hat+exp(x_j'*theta_nl)^2*x_j*x_j'/n;
    B_hat=B_hat+u_hat_j^2*exp(x_j'*theta_nl)^2*x_j*x_j'/n;
end
V_hat=inv(A_hat)*B_hat*inv(A_hat);
