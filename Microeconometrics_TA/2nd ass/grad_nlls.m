function value=grad_nlls(theta)

global n k x y

sum_grad=zeros(k,1);

for i=1:n
    sum_grad=sum_grad+(y(i,1)-exp(x(i,:)*theta))*exp(x(i,:)*theta)*x(i,:)'/((-1)*n);
end

value=sum_grad;

end
