%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%SIMULATED MAXIMUM LIKELIHOOD(SML)%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BHHH ITERATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
load('DGP_nonzero_corr');
[N,J]=size(Big_Y(:,:,1));R=length(Big_mu(:,:,1));
true_beta_2=[0.5;1;1];true_beta_3=[0.5;1;1];c1=1.33;c2=1;
true_theta=[true_beta_2;true_beta_3;c1;c2];K=length(true_theta);
MC_N=500;
lambda=0.5;
theta_SML=zeros(K,MC_N);
R=30;Big_mu_R=Big_mu(1:R,:,:);
for m=1:MC_N
    X=Big_X(:,:,m);mu=Big_mu_R(:,:,m);Y=Big_Y(:,:,m);
    iter=1;
    initial_theta_SML=true_theta+lambda*BHHH(Y,X,mu,tilde_L,true_theta);
    theta_step_SML=initial_theta_SML;
    criterion=norm(initial_theta_SML-true_theta);
    while criterion>0.0001
        theta_next_SML=theta_step_SML+lambda*BHHH(Y,X,mu,tilde_L,theta_step_SML);
        iter=iter+1
        criterion=norm(theta_next_SML-theta_step_SML)
        theta_step_SML=theta_next_SML;
        if iter>50
            break
        end
      
    end
    
    theta_SML(:,m)=theta_next_SML;

end


