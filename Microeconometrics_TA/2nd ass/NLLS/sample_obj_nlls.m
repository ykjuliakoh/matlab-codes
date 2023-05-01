function value=sample_obj_nlls(theta)


global N x y

value0=0;

for j=1:N
    x_j=x(j,:)';
    y_j=y(j,1);
    u_j=(y_j-exp(x_j(1)*theta(1)+x_j(2)*theta(2)+x_j(3)*theta(3)) );
    value0=value0+u_j^2;
end

value=value0;