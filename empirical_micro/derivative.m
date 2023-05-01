function derivative_matrix=derivative(X,L,mu,j,theta_step)

[N,J]=size(X);R=length(mu);
K=length(theta_step);
beta_2=theta_step(1:3,:);beta_3=theta_step(4:6,:);
X_2=[ones(N,1),X(:,1),X(:,2)];
X_3=[ones(N,1),X(:,1),X(:,3)];
tilde_V_21=X_2*beta_2;tilde_V_31=X_3*beta_3;
tilde_V_12=-tilde_V_21;tilde_V_32=tilde_V_31-tilde_V_21;
tilde_V_13=-tilde_V_31;tilde_V_23=-tilde_V_32;

tilde_V_full=zeros(N,J-1,J);
tilde_V_full(:,:,1)=[tilde_V_21,tilde_V_31];
tilde_V_full(:,:,2)=[tilde_V_12,tilde_V_32];
tilde_V_full(:,:,3)=[tilde_V_13,tilde_V_23];

L_j=L(:,:,j);
c_aa=L_j(1,1);c_ab=L_j(2,1);c_bb=L_j(2,2);

matrix=zeros(N,K);
if j==1
    for i=1:N
        A=-tilde_V_full(i,1,j)/c_aa;
        eta_r=norminv(normcdf(A)*mu(:,j));
        B=-(tilde_V_full(i,2,j)*ones(R,1)+c_ab*eta_r)/c_bb;
        Dev_1=(-normpdf(A)/c_aa)*normcdf(B)+(normcdf(A)*normpdf(A)*(-c_ab/c_bb)*(-1/c_aa))*(normpdf(B).*(normpdf(eta_r).^-1).*mu(:,j));
        Dev_2=(-normcdf(A)/c_bb)*normpdf(B);
        matrix(i,1:3)=(sum(Dev_1)/R)*X_2(i,:);
        matrix(i,4:6)=(sum(Dev_2)/R)*X_3(i,:);
    end
elseif j==2
    for i=1:N
        
        A=-tilde_V_full(i,1,j)/c_aa;
        eta_r=norminv(normcdf(A)*mu(:,j));
        B=-(tilde_V_full(i,2,j)*ones(R,1)+c_ab*eta_r)/c_bb;
        Dev_1=(-normpdf(A)/c_aa)*normcdf(B)+normcdf(A)*normpdf(B).*((-1/c_bb)*ones(R,1)+((-c_ab/c_bb)*(-1/c_aa)*normpdf(A))*((normpdf(eta_r).^-1).*mu(:,j)));
        Dev_2=(-normcdf(A)/c_bb)*normpdf(B);
        matrix(i,1:3)=(sum(Dev_1)/R)*(-X_2(i,:));
        matrix(i,4:6)=(sum(Dev_2)/R)*X_3(i,:);
    end
    
else
    for i=1:N
        A=-tilde_V_full(i,1,j)/c_aa;
        eta_r=norminv(normcdf(A)*mu(:,j));
        B=-(tilde_V_full(i,2,j)*ones(R,1)+c_ab*eta_r)/c_bb;
        Dev_1=(-normcdf(A)/c_bb)*normpdf(B);
        Dev_2=(-normpdf(A)/c_aa)*normcdf(B)+normcdf(A)*normpdf(B).*((-1/c_bb)*ones(R,1)+((-c_ab/c_bb)*normpdf(A)*(-1/c_aa))*((normpdf(eta_r).^-1).*mu(:,j)));
        matrix(i,1:3)=(sum(Dev_1)/R)*X_2(i,:);
        matrix(i,4:6)=(sum(Dev_2)/R)*(-X_3(i,:));
    end
end
derivative_matrix=matrix;
end

        