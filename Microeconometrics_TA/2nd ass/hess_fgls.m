function value=hess_fgls(theta)

global k_temp uhat_sq x_temp n_temp

sum_hess=zeros(k_temp,k_temp);

for i=1:n_temp
    sum_hess=sum_hess+(uhat_sq(i,1)-2*exp(x_temp(i,:)*theta))*exp(x_temp(i,:)*theta)*x_temp(i,:)'*x_temp(i,:)/((-1)*n_temp);
end

value=sum_hess;

end