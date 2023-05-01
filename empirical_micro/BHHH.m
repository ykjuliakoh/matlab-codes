function vector_valued=BHHH(Y,X,mu,L,theta)

[N,J]=size(X);K=length(theta);
theta_no_c=theta(1:6,1);
simulated_p=GHK_simulator(X,mu,L,theta_no_c);%(N,J)

derivative_matrix=full_derivative(X,mu,L,theta);%(K,NJ)
S_N=zeros(K,N);
for i=1:N
    Y_i=Y(i,:);
    simulated_p_i=simulated_p(i,:);
    derivative_i=derivative_matrix(:,3*i-2:3*i);
    S_N(:,i)=derivative_i*(Y_i./simulated_p_i)';
end
initial_hessian_sum=zeros(K,K);
for i=1:N
    S_i=S_N(:,i);
    hessian_sum=initial_hessian_sum+S_i*S_i';
    initial_hessian_sum=hessian_sum;
end
gradient=sum(S_N,2);
hessian=hessian_sum;
vector_valued=inv(hessian)*gradient;



