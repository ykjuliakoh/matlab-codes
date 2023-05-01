function value=hess_nlls(theta)

global n k x y

sum_hess=zeros(k,k);

for i=1:n
    sum_hess=sum_hess+(y(i,1)-2*exp(x(i,:)*theta))*exp(x(i,:)*theta)*x(i,:)'*x(i,:)/((-1)*n);
end

value=sum_hess;

end