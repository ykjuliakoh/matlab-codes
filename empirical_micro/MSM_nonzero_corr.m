clear;
load('DGP_nonzero_corr');
[N,J]=size(Big_Y(:,:,1));R=length(Big_mu(:,:,1));
true_beta_2=[0.5;1;1];true_beta_3=[0.5;1;1];c1=1.33;c2=1;
true_theta=[true_beta_2;true_beta_3;c1;c2];
K=length(true_theta);
MC_N=500;
lambda=0.5;
R=30;Big_mu_R=Big_mu(1:R,:,:);
theta_MSM_R50=zeros(K,MC_N);
for m=1:MC_N
    X=Big_X(:,:,m);mu=Big_mu_R(:,:,m);Y=Big_Y(:,:,m);
    iter=1;
    initial_theta_MSM=true_theta+lambda*GN(Y,X,mu,tilde_L,true_theta);
    theta_step_MSM=initial_theta_MSM;
    criterion=norm(initial_theta_MSM-true_theta);
    while criterion>0.00001
        theta_next_MSM=theta_step_MSM+lambda*GN(Y,X,mu,tilde_L,theta_step_MSM);
        iter=iter+1;
        criterion=norm(theta_next_MSM-theta_step_MSM);
        theta_step_MSM=theta_next_MSM;
        if iter>2000
            break
        end

    end
    
    theta_MSM_R50(:,m)=theta_next_MSM;

end

