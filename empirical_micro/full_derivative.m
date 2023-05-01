
function derivative_matrix=full_derivative(X,mu,L,theta)
[N,J]=size(X);K=length(theta);
X_2=[ones(N,1),X(:,1),X(:,2)];
X_3=[ones(N,1),X(:,1),X(:,3)];
beta_2=theta(1:3,1);beta_3=theta(4:6,1);
c1=theta(7,1);c2=theta(8,1);

partial_beta=zeros(K,N*J);
for j=1:J
    L_j=L(:,:,j);
    c_aa=L_j(1,1);c_ab=L_j(2,1);c_bb=L_j(2,2);
    if j==1
        tilde_V_first=X_2*beta_2;tilde_V_second=X_3*beta_3;
        for i=1:N
            X_2i=X_2(i,:);
            X_3i=X_3(i,:);
            mu_i=mu(:,3*i-2:3*i);mu_ij=mu_i(:,j);R=length(mu_ij);
            tilde_V_first_i=tilde_V_first(i,1);
            tilde_V_second_i=tilde_V_second(i,1);
            a=-tilde_V_first_i/c_aa;
            eta_r=norminv(mu_ij*normcdf(a));
            b=-(tilde_V_second_i*ones(R,1)+c_ab*eta_r)./c_bb;
            partial_beta(1:3,3*i-2)=(normpdf(a)*(-1/c_aa))*sum(normcdf(a)*normpdf(b)*(-c_ab/c_bb).*(normpdf(eta_r).^-1).*mu_ij+normcdf(b))/R*X_2i';
            partial_beta(4:6,3*i-2)=(normcdf(a)*(-1/c_bb))*sum(normpdf(b))/R*X_3i';
            partial_beta(7,3*i-2)=(normcdf(a)*(-1/c_bb))*sum(normpdf(b).*eta_r)/R;
            partial_beta(8,3*i-2)=(normcdf(a)/c_bb)*sum(normpdf(b).*(-b))/R;
        end
    elseif j==2
        tilde_V_first=-X_2*beta_2;tilde_V_second=X_3*beta_3-X_2*beta_2;
        for i=1:N
            X_2i=X_2(i,:);
            X_3i=X_3(i,:);
            mu_i=mu(:,3*i-2:3*i);mu_ij=mu_i(:,j);R=length(mu_ij);
            tilde_V_first_i=tilde_V_first(i,1);
            tilde_V_second_i=tilde_V_second(i,1);
            a=-tilde_V_first_i/c_aa;
            eta_r=norminv(mu_ij*normcdf(a));
            b=-(tilde_V_second_i*ones(R,1)+c_ab*eta_r)./c_bb;
            partial_beta(1:3,3*i-1)=-((normpdf(a)*(-1/c_aa)*sum(normcdf(b)+normcdf(a)*normpdf(b)*(-c_ab/c_bb).*(normpdf(eta_r).^-1).*mu_ij)/R)+(normcdf(a)*(-1/c_bb)*sum(normpdf(b))/R))*X_2i';
            partial_beta(4:6,3*i-1)=(normcdf(a)*(-1/c_bb)*sum(normpdf(b))/R)*X_3i';
            partial_beta(7,3*i-1)=(normcdf(a)*(1/c_bb))*sum(normpdf(b).*eta_r)/R;
            partial_beta(8,3*i-1)=(normcdf(a)*(1/c_bb))*sum(normpdf(b).*(-b))/R;
        end
    else 
        tilde_V_first=-X_3*beta_3;
        tilde_V_second=X_2*beta_2-X_3*beta_3;
        for i=1:N
            X_2i=X_2(i,:);
            X_3i=X_3(i,:);
            mu_i=mu(:,3*i-2:3*i);mu_ij=mu_i(:,j);R=length(mu_ij);
            tilde_V_first_i=tilde_V_first(i,1);
            tilde_V_second_i=tilde_V_second(i,1);
            a=-tilde_V_first_i/c_aa;
            eta_r=norminv(mu_ij*normcdf(a));
            b=-(tilde_V_second_i*ones(R,1)+c_ab*eta_r)./c_bb;
            partial_beta(1:3,3*i)=normcdf(a)*(-1/c_bb)*sum(normpdf(b))/R*X_2i';
            partial_beta(4:6,3*i)=-((normpdf(a)*(-1/c_aa)*sum(normcdf(b)+normcdf(a)*normpdf(b)*(-c_ab/c_bb).*(normpdf(eta_r).^-1).*mu_ij)/R)+(-1/c_bb)*normcdf(a)*sum(normpdf(b))/R)*X_3i';
            partial_beta(7,3*i)=((c1*(c1^2+c2^2)^(-1/2)*normpdf(a)*(-a/c_aa))*sum(normcdf(b))/R)+(((c1^2+c2^2)^(-3/2)*(c1^3+c1*c2^2-c2^2)*normcdf(a)*(-1/c_bb))*sum(normpdf(b).*eta_r)/R)+((-c1*c2*(c1^2+c2^2)^(-3/2)*normcdf(a)*(1/c_bb))*sum(normpdf(b).*(-b))/R);
            partial_beta(8,3*i)=((c2*(c1^2+c2^2)^(-1/2)*normpdf(a)*(-a/c_aa))*sum(normcdf(b))/R)+(((c1^2+c2^2)^(-3/2)*(c2^3+c1^2*c2+c1*c2)*normcdf(a)*(-1/c_bb))*sum(normpdf(b).*eta_r)/R)+((c1^2*(c1^2+c2^2)^(-3/2)*normcdf(a)*(1/c_bb))*sum(normpdf(b).*(-b))/R);
        end
    end
end

derivative_matrix=partial_beta;
end

        
            