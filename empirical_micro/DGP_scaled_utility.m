%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%DATA GENERATING PROCESS%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
rng(1);
MC_N=500;%NUMBER OF MONTECARLO SIMULATION
N=500;J=3;
Big_X=normrnd(0,1,[N,3,MC_N]);
MU=zeros(J-1,1);
%SCALED UTILITY%normalized the utility for an alternative 1
c1=1.33;c2=1;
beta_2=[0.5,1,1]';beta_3=[0.5,1,1]';
Big_U=zeros(N,J,MC_N);Big_Y=zeros(N,J,MC_N);
R=100;%replication for GHK simulator

Big_mu=rand(R,N*J,MC_N);

for m=1:MC_N
    rng(m);
    epsilon_m=mvnrnd(MU,eye(J-1),N);
    X_m=Big_X(:,:,m);
    X_2=[ones(N,1),X_m(:,1),X_m(:,2)];
    X_3=[ones(N,1),X_m(:,1),X_m(:,3)];
    V_2=X_2*beta_2;V_3=X_3*beta_3;
    u_2=V_2+epsilon_m(:,1);
    u_3=V_3+c1*epsilon_m(:,1)+c2*epsilon_m(:,2);   
    Big_U(:,2,m)=u_2;Big_U(:,3,m)=u_3;
    for i=1:N
        Big_U_i=Big_U(i,:,m);
        [M,I]=max(Big_U_i);
        Big_Y(i,I,m)=1;
    end

end
sigma_tilde_1=zeros(J-1,J-1);
sigma_tilde_1(1,1)=1;sigma_tilde_1(2,1)=c1;
sigma_tilde_1(1,2)=c1;sigma_tilde_1(2,2)=c1^2+c2^2;
L1=chol(sigma_tilde_1,'lower');
L=zeros(J,J);
L(2:J,2:J)=L1;sigma=L*L';
identity_matrix=eye(J-1);
M_1=[-ones(J-1,1),identity_matrix];M_3=[identity_matrix,-ones(J-1,1)];
M_2=[identity_matrix(:,1),-ones(J-1,1),identity_matrix(:,2)];
M=zeros(J-1,J,J);M(:,:,1)=M_1;M(:,:,2)=M_2;M(:,:,3)=M_3;
tilde_sigma=zeros(J-1,J-1,J);tilde_L=zeros(J-1,J-1,J);
for j=1:J
    tilde_sigma(:,:,j)=M(:,:,j)*sigma*M(:,:,j)';
    tilde_L(:,:,j)=chol(tilde_sigma(:,:,j),'lower');
end

%%%%%SIMULATION FOR GII%%%%%%
M1=10;M2=300;
BIG_EPSILON=normrnd(0,1,[N*2,M1+M2,MC_N]);
Big_Z=zeros(N,11,MC_N);
for m=1:MC_N
    X=Big_X(:,:,m);
    x1=X(:,1);x2=X(:,2);x3=X(:,3);
    Big_Z(:,:,m)=[ones(N,1),x1,x2,x3,x1.^2,x2.^2,x3.^2,x1.*x2,x1.*x3,x2.*x3,x1.*x2.*x3];
end

Big_theta_hat=zeros(22,MC_N);
for m=1:MC_N
    Y=Big_Y(:,:,m);Z=Big_Z(:,:,m);
    Big_theta_hat(1:11,m)=inv(Z'*Z)*Z'*Y(:,2);
    Big_theta_hat(12:22,m)=inv(Z'*Z)*Z'*Y(:,3);
end


save('DGP_nonzero_corr','Big_X','Big_Y','Big_theta_hat','BIG_EPSILON','Big_Z','tilde_L','Big_mu')