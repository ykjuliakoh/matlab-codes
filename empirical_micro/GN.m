function vector_valued=GN(Y,X,mu,L,theta)
[N,J]=size(X);
theta_no_c=theta(1:6,1);
simulated_p_theta=GHK_simulator(X,mu,L,theta_no_c);
simulated_p=zeros(N*J,1);
Y_NJ=zeros(N*J,1);
for i=1:N
    simulated_p(3*i-2:3*i,1)=simulated_p_theta(i,:)';
    Y_NJ(3*i-2:3*i,1)=Y(i,:)';
end

derivative_matrix=full_derivative(X,mu,L,theta);

W=(derivative_matrix*(simulated_p.^-1))'*(derivative_matrix*(simulated_p.^-1));
gradient=-derivative_matrix*W*(Y_NJ-simulated_p);
hessian=derivative_matrix*W*derivative_matrix';

vector_valued=inv(-hessian)*gradient;