% using function "grad_probit.m" and "hess_probit.m"
theta_ini=[1,1]';
% Or try anther initial value

%test_g=grad_probit;
%test_h=hess_probit;

theta_old=theta_ini;
theta_new=theta_old-inv(hess_probit(theta_old)) * grad_probit(theta_old);

iter=1;
dist= norm(theta_new-theta_old);

while dist>0.00001
  theta_old=theta_new -inv(hess_probit(theta_new))*grad_probit(theta_new);
  iter=iter+1;
  dist= norm(theta_new-theta_old);
  theta_new=theta_old;
   if iter>2000
       break
   end
end
theta_nl=theta_new;