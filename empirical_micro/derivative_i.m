function matrix=derivative_i(X,L,mu,theta)
K=length(theta);J=length(X);
X_2i=[1,X(:,1),X(:,2)];X_3i=[1,X(:,1),X(:,3)];
beta2=theta(1:3,1);beta3=theta(4:6,1);

matrix_j=zeros(K,J);

for j=1:J
    L_j=L(:,:,j);mu_j=mu(:,j);R=length(mu_j);
    c_aa=L_j(1,1);c_ab=L_j(2,1);c_bb=L_j(2,2);
    if j==1
        V_tilde_first_i=X_2i*beta2;
        V_tilde_second_i=X_3i*beta3;
        a=-V_tilde_first_i/c_aa;
        eta_r=norminv(mu_j*normcdf(a));
        b=-(V_tilde_second_i*ones(R,1)+c_ab*eta_r)./c_bb;
        matrix_j(1:3,j)=X_2i'*(normpdf(a)*(-1/c_aa))*(sum((normcdf(a)*normpdf(b)*(-c_ab/c_bb).*(normpdf(eta_r).^-1).*mu_j)+normcdf(b))/R);
        matrix_j(4:6,j)=X_3i'*(normcdf(a)*(-1/c_bb))*(sum(normpdf(b))/R);
    elseif j==2
        V_tilde_first_i=-X_2i*beta2;
        V_tilde_second_i=X_3i*beta3-X_2i*beta2;
        a=-V_tilde_first_i/c_aa;
        eta_r=norminv(mu_j*normcdf(a));
        b=-(V_tilde_second_i*ones(R,1)+c_ab*eta_r)./c_bb;
        matrix_j(1:3,j)=-X_2i'*((normpdf(a)*(-1/c_aa))*(sum(normcdf(b)+normcdf(a)*normpdf(b)*(-c_ab/c_bb).*(normpdf(eta_r).^-1).*mu_j)/R)+normcdf(a)*(-1/c_bb)*(sum(normpdf(b))/R));
        matrix_j(4:6,j)=X_3i'*normcdf(a)*(-1/c_bb)*(sum(normpdf(b))/R);
    else
        V_tilde_first_i=-X_3i*beta3;
        V_tilde_second_i=X_2i*beta2-X_3i*beta3;
        a=-V_tilde_first_i/c_aa;
        eta_r=norminv(mu_j*normcdf(a));
        b=-(V_tilde_second_i*ones(R,1)+c_ab*eta_r)./c_bb;
        matrix_j(1:3,j)=X_2i'*normcdf(a)*(-1/c_bb)*(sum(normpdf(b))/R);
        matrix_j(4:6,j)=-X_3i'*(normpdf(a)*(-1/c_aa)*(sum(normcdf(b)+normcdf(a)*normpdf(b)*(-c_ab/c_bb).*(normpdf(eta_r).^-1).*mu_j)/R)+(-1/c_bb)*normcdf(a)*(sum(normpdf(b))/R));
    end
end

matrix=matrix_j;

        
        