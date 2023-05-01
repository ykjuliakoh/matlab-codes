function matrix=hess_nlls(theta)

% theta is a column vector

global N x y K

matrix_sum=zeros(K,K);

for j=1:N
    x_j=x(j,:)';
    y_j=y(j,1);
    hess_j= -(y_j-2*exp(x_j'*theta))*exp(x_j'*theta)*x_j*x_j';
    matrix_sum=matrix_sum + hess_j;
end

matrix=matrix_sum/N;