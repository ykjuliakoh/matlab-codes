function matrix=hess_probit(theta)

% theta is a column vector

global N x y K

matrix_sum=zeros(K,K);

for j=1:N
    x_j=x(j,:);
    y_j=y(j,1);
    hess_j=normcdf((x_j)*theta)^2*((y_j)-normcdf((x_j)*theta))^2/(normcdf((x_j)*theta)^2*(1-normcdf(x_j*theta))^2)*(x_j')*(x_j);
    matrix_sum=matrix_sum + hess_j;
end

matrix=matrix_sum/N;