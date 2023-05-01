function vector=gradient_GN(theta)
global K N simulated_p Big_y
%derivative matrix is (J,K) matrix
initial_gradient_sum=zeros(K,1);W=instrument_matrix(theta);
for i=1:N
    derivative_matrix_i=derivative(i,theta);
    simulated_p_i=simulated_p(i,:)';y_i=Big_y(i,:)';
    moment_condition_i=y_i-simulated_p_i;
    W_i=W(i,:);%(1,K)
    gradient_i=derivative_matrix_i'*(W_i*W_i')*moment_condition_i;
    gradient_sum=initial_gradient_sum+gradient_i;
    initial_gradient_sum=gradient_sum;    
end
vector=-gradient_sum;
end


