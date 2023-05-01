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


