function value=grad_fgls(theta)

global n_temp k_temp uhat_sq x_temp

sum_grad=zeros(k_temp,1);

for i=1:n_temp
    sum_grad=sum_grad+(uhat_sq(i,1)-exp(x_temp(i,:)*theta))*exp(x_temp(i,:)*theta)*x_temp(i,:)'/((-1)*n_temp);
end

value=sum_grad;

end
