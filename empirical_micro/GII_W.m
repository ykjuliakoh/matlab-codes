clear;
load('GII_simul500');
MC_N=500;%MONTE CARLO REPLICATION NUMBER
true_beta_2=[0.5;1;1];true_beta_3=[0.5;1;1];
true_beta=[true_beta_2;true_beta_3];K=length(true_beta);
alpha=.5;

results_GII_W_QN_step1=zeros(K,MC_N);
for m=1:MC_N
    X=Big_X(:,:,m);theta_hat=Big_theta_hat(:,m);
    Y=Big_Y(:,:,m);
    Z=Big_Z(:,:,m);simul_e_step1=BIG_EPSILON(:,1:10,m);
    initial_results=GII_WALD(X,Y,Z,theta_hat,simul_e_step1,true_beta);
    initial_gradient=initial_results(:,1);
    initial_hessian=initial_results(:,2:7);
    initial_beta=true_beta+alpha*inv(-initial_hessian)*initial_gradient;
    beta_W=initial_beta;
    criterion=norm(initial_beta-true_beta);
    iter=1;

    while criterion>0.00001
        results=GII_WALD_QN(X,Y,Z,theta_hat,simul_e_step1,beta_W,initial_gradient,initial_hessian);
        iter=iter+1
        gradient=results(:,1);hessian=results(:,2:7);
        beta_next_W=beta_W+alpha*inv(-hessian)*gradient;
        criterion=norm(beta_next_W-beta_W)
        beta_W=beta_next_W;
        initial_gradient=gradient;initial_hessian=hessian;
        if iter>2000
            break
        end

    end
    results_GII_W_QN_step1(:,m)=beta_W;
end
