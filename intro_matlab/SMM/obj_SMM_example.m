function value=obj_SMM_example(b)

global x y N S v_tilde_SN

beta=[b(1), b(2)]';
value0=0;
for i=1:N
    x_i=x(i,:)';

    v_tilde_SN_i=v_tilde_SN(:,i);  %  y_star_tilde_i_S is S*1 vector      
    p_tilde_i=sum(v_tilde_SN_i>-1*x_i'*beta)/S;
    w_i=x_i*normpdf(-x_i'*beta)/((1-p_tilde_i)*p_tilde_i);
    value0=value0+(y(i,1)-p_tilde_i)*w_i;
end

value=(value0/N)'*(value0/N);