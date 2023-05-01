clear;
load('GII_simul500');
MC_N=500;%MONTE CARLO REPLICATION NUMBER
true_beta_2=[0.5;1;1];true_beta_3=[0.5;1;1];
true_beta=[true_beta_2;true_beta_3];K=length(true_beta);
alpha=.5;

results_GII_LR_QN_step1=zeros(K,MC_N);
for m=1:MC_N
    X=Big_X(:,:,m);theta_hat=Big_theta_hat(:,m);
    Y=Big_Y(:,:,m);
    Z=Big_Z(:,:,m);epsilon=BIG_EPSILON(:,1:10,m);
    initial_results=GII_LR(X,Y,Z,theta_hat,epsilon,true_beta);
    initial_gradient=initial_results(:,1);
    initial_hessian=initial_results(:,2:7);
    initial_beta=true_beta+alpha*inv(-initial_hessian)*initial_gradient;
    beta_LR=initial_beta;
    criterion=norm(initial_beta-true_beta);
    iter=1;

    while criterion>0.00001
        results=GII_LR_QN(X,Y,Z,theta_hat,epsilon,beta_LR,initial_gradient,initial_hessian);
        iter=iter+1
        gradient=results(:,1);hessian=results(:,2:7);
        beta_next_LR=beta_LR+alpha*inv(-hessian)*gradient;
        criterion=norm(beta_next_LR-beta_LR)
        beta_LR=beta_next_LR;
        initial_gradient=gradient;initial_hessian=hessian;
        if iter>2000
            break
        end

    end
    results_GII_LR_QN_step1(:,m)=beta_LR;
end
%%%%%%%%%%%SECOND STEP NR%%%%%%%%%%%%%%%%%%%
clear;
load('GII_simul500');
load('results_GII_LR_QN1.mat');
GII_LR=results_GII_LR_QN_step1(:,~all(isnan(results_GII_LR_QN_step1)));
[K,M_LR]=size(GII_LR);
beta_GII_LR2=zeros(6,M_LR);
for m=1:M_LR
    X=Big_X(:,:,m);theta_hat=Big_theta_hat(:,m);
    Y=Big_Y(:,:,m);
    Z=Big_Z(:,:,m);epsilon=BIG_EPSILON(:,11:110,m);
    beta_0=GII_LR(:,m);
    
    beta_GII_LR2(:,m)=beta_0+GII_LR2(X,Y,Z,theta_hat,epsilon,beta_0);
    
end

results_2=beta_GII_LR2(:,~all(isnan(beta_GII_LR2)));
[K,M_LR2]=size(results_2);
