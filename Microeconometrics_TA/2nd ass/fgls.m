global x k_temp uhat_sq x_temp n_temp

%%%% Run FGLS %%%%%%%%%%%%%%%%%%%%%
%%% Obtain theta_ols and uhat %%%
x_temp=[x(:,2),x(:,3)];
[n_temp,k_temp]=size(x_temp);

sum_xx=zeros(k_temp,k_temp);
sum_xy=zeros(k_temp,1);

for i=1:n
    sum_xx=sum_xx+x_temp(i,:)'*x_temp(i,:)/n;
    sum_xy=sum_xy+x_temp(i,:)'*y(i,1)/n;
end

theta_ols=inv(sum_xx)*sum_xy;

uhat_=y-x_temp*theta_ols;
uhat_sq=uhat.^2;

%%% Obtain alpha_hat: Regress uhat_sq on the x nonlinearlly %%%

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
