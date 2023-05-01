function vector=grad_nlls(theta)

% theta is a column vector

global N x y K

vector_sum=zeros(K,1);

for j=1:N
    x_j=x(j,:)';
    y_j=y(j,1);
    derivative_j=(y_j-exp(x_j'*theta))*(-1)*exp(x_j'*theta)*x_j;
    vector_sum=vector_sum+derivative_j;
end

vector=vector_sum/N;