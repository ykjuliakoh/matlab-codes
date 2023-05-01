%%% Consolidated coding %%%

%%% OLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,k]=size(x);

sum_xx=zeros(k,k);
sum_xy=zeros(k,1);

for i=1:n
    sum_xx=sum_xx+x(i,:)'*x(i,:)/n;
    sum_xy=sum_xy+x(i,:)'*y(i,1)/n;
end

theta_ols=inv(sum_xx)*sum_xy;

uhat=y-x*theta_ols;

%%% GLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_uxx=zeros(k,k);

for i=1:n
    sum_uxx=sum_uxx+uhat(i,1)^2*x(i,:)'*x(i,:)/n;
end

av_ols=inv(sum_xx)*sum_uxx*inv(sum_xx);
ste_ols=sqrt(diag(av_ols)./n);

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

%%%% FGLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global x k_temp uhat_sq x_temp n_temp

x=[ones,x1,x2];

[n,k]=size(x);

%%% step1: Obtain theta_ols and uhat %%%
x_temp=[x1,x2];
[n_temp,k_temp]=size(x_temp);

sum_xx=zeros(k_temp,k_temp);
sum_xy=zeros(k_temp,1);

for i=1:n
    sum_xx=sum_xx+x_temp(i,:)'*x_temp(i,:)/n;
    sum_xy=sum_xy+x_temp(i,:)'*y(i,1)/n;
end

theta_ols=inv(sum_xx)*sum_xy;

uhat=y-x_temp*theta_ols;
uhat_sq=uhat.^2;

%%% step2: Obtain alpha_hat: Regress uhat_sq on the x "NONLINEARLLY" %%%

iter=1;

theta_old=theta_ols;
theta_new=theta_old-inv(hess_fgls(theta_old))*grad_fgls(theta_old);
dist=norm(theta_new-theta_old);

theta_old=theta_new;

while dist>0.0000001
    theta_new=theta_old-inv(hess_fgls(theta_old))*grad_fgls(theta_old);
    iter=iter+1;
    dist=norm(theta_new-theta_old);
    theta_old=theta_new;
    
    if iter>2000
        break
    end
    
end

alpha_hat=theta_new;

%%% Create predictions sigma_sq_hat %%%
sigma_sq_hat=exp(x_temp*alpha_hat);

%%% Obtain theta_fgls by substituting sigma_sq with sigma_sq_hat %%%
sum_xx_fgls=zeros(k,k);
sum_xy_fgls=zeros(k,1);

for i=1:n
    sum_xx_fgls=sum_xx_fgls+x(i,:)'*x(i,:)/(sigma_sq_hat(i,1)*n);
    sum_xy_fgls=sum_xy_fgls+x(i,:)'*y(i,1)/(sigma_sq_hat(i,1)*n);
end

theta_fgls=sum_xx_fgls*sum_xy_fgls;

%%% Obtain av_fgls and ste_fgls %%%
uhat_fgls=y-x*theta_fgls;

sum_uxx_fgls=zeros(k,k);

for i=1:n
    sum_uxx_fgls=sum_uxx_fgls+uhat_fgls(i,1)^2*x(i,:)'*x(i,:)/(sigma_sq_hat(i,1)*n);
end

av_fgls=inv(sum_xx_fgls)*sum_uxx_fgls*inv(sum_xx_fgls);
ste_fgls=sqrt(diag(av_fgls)./n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CIV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
ste_civ=sqrt(diag(av_civ)./n);

%%% TSLS: CIV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Note: Use AV_civ, not SE of the 2nd stage!

%%% GIV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%% TSLS: GIV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Note: Use AV_GIV, not SE of the 2nd stage!

%%% NLLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,k]=size(x);

global n k x y

iter=1;
theta_ols=inv(x'*x)*x'*y;

theta_old=theta_ols;
theta_new=theta_old-inv(hess_nlls(theta_old))*grad_nlls(theta_old);
dist=norm(theta_new-theta_old);

theta_old=theta_new;

while dist>0.0000001
    theta_new=theta_old-inv(hess_nlls(theta_old))*grad_nlls(theta_old);
    iter=iter+1;
    dist=norm(theta_new-theta_old);
    theta_old=theta_new;
    
    if iter>2000
        break
    end
    
end

theta_nlls=theta_new;

uhat=y-x*theta_nlls;

global n k

uhat=y-x*theta_nlls;

%%% hetero %%%
sum_xx=zeros(k,k);
sum_uxx=zeros(k,k);

for i=1:n
    sum_xx=sum_xx+exp(2*x(i,:)*theta_nlls)*x(i,:)'*x(i,:)/n;
    sum_uxx=sum_uxx+uhat(i,1)^2*exp(2*x(i,:)*theta_nlls)*x(i,:)'*x(i,:)/n;
end

Ahat=sum_xx;
Bhat=sum_uxx;

av_nlls_hetero=inv(Ahat)*Bhat*inv(Ahat);
ste_nlls_hetero=sqrt(diag(av_nlls_hetero)./n);


%%% homo %%%
sigma_sq=zeros(1,1);

for i=1:n
    sigma_sq=sigma_sq+uhat(i,1)^2/n;
end

av_nlls_homo=sigma_sq*inv(Ahat);
ste_nlls_homo=sqrt(diag(av_nlls_homo)./n);

%%% General t-test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global n

q_grad=[0,1,1,1]';

%%% hetero %%%
v_q_hetero=q_grad'*av_nlls_hetero*q_grad;

t_q_hetero=sqrt(n)*(theta_nlls(2)+theta_nlls(3)+theta_nlls(4))/v_q_hetero;

%%% homo %%%
v_q_homo=q_grad'*av_nlls_homo*q_grad;

t_q_homo=sqrt(n)*(theta_nlls(2)+theta_nlls(3)+theta_nlls(4))/v_q_homo;

%%% Stata NLLS command: nl(reg model) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, under the homoskedasticity: nl(y=exp({ones}+{x1}*x1+{x2}*x2+{x3}*x3))
% Here, under the heteroskedasticity:
% nl(y=exp({ones}+{x1}*x1+{x2}*x2+{x3}*x3)), vce(r)

