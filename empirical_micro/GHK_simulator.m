
function matrix=GHK_simulator(X,mu,L,theta)
[N,J]=size(X);
X_2=[ones(N,1),X(:,1),X(:,2)];X_3=[ones(N,1),X(:,1),X(:,3)];
beta2=theta(1:3,1);beta3=theta(4:6,1);

probability_matrix=zeros(N,J);
for j=1:J
    L_j=L(:,:,j);
    c_aa=L_j(1,1);c_ab=L_j(2,1);c_bb=L_j(2,2);
    if j==1
        tilde_V_first=X_2*beta2;
        tilde_V_second=X_3*beta3;
    elseif j==2
        tilde_V_first=-X_2*beta2;
        tilde_V_second=X_3*beta3-X_2*beta2;
    else
        tilde_V_first=-X_3*beta3;
        tilde_V_second=X_2*beta2-X_3*beta3;
    end
    for i=1:N
        tilde_V_first_i=tilde_V_first(i,1);
        tilde_V_second_i=tilde_V_second(i,1);
        k=normcdf(-tilde_V_first_i/c_aa);
        mu_i=mu(:,3*i-2:3*i);mu_ij=mu_i(:,j);R=length(mu_ij);
        eta_r=norminv(mu_ij*k);
        g_r=normcdf(-(tilde_V_second_i*ones(R,1)+c_ab*eta_r)./c_bb);
        probability_matrix(i,j)=sum(k*g_r)/R;
    end
end
matrix=probability_matrix;
end

    
