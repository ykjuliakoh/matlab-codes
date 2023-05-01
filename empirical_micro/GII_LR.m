function gradient_hessian=GII_LR(X,Y,Z,theta_hat,epsilon,beta)
lambda_1=0.03;
[N,J]=size(X);M=10;p=length(theta_hat);K=length(beta);
epsilon_2=epsilon(1:N,:);epsilon_3=epsilon(N+1:2*N,:);
X_2=[ones(N,1),X(:,1),X(:,2)];X_3=[ones(N,1),X(:,1),X(:,3)];
beta_2=beta(1:3,1);beta_3=beta(4:6,1);
stacked_Y=[Y(:,2);Y(:,3)];
BZ=[Z,zeros(N,11);zeros(N,11),Z];
multiplied_BZ=inv(BZ'*BZ)*BZ';
simulated_u2=zeros(N,M);simulated_u3=zeros(N,M);
for m=1:M
    simulated_u2(:,m)=X_2*beta_2+epsilon_2(:,m);
    simulated_u3(:,m)=X_3*beta_3+epsilon_3(:,m);
end
smoothed_u2=zeros(N,M);smoothed_u3=zeros(N,M);
for m=1:M
    smoothed_u2(:,m)=exp(simulated_u2(:,m)/lambda_1)./(ones(N,1)+exp(simulated_u2(:,m)/lambda_1)+exp(simulated_u3(:,m)/lambda_1));
    smoothed_u3(:,m)=exp(simulated_u3(:,m)/lambda_1)./(ones(N,1)+exp(simulated_u2(:,m)/lambda_1)+exp(simulated_u3(:,m)/lambda_1));
end
theta_simul=multiplied_BZ*[smoothed_u2;smoothed_u3];
theta_bar=sum(theta_simul,2)/M;

partial_g2=zeros(N,K);partial_g3=zeros(N,K);
for k=1:K
    stacked_X=[ones(N,1),X(:,1),X(:,2),ones(N,1),X(:,1),X(:,3)];
    denominator=(ones(N,M)+exp(simulated_u2./lambda_1)+exp(simulated_u3./lambda_1)).^2;
    coeff_A=(1/lambda_1)*(ones(N,M)+exp(simulated_u3./lambda_1)).*exp(simulated_u2./lambda_1);
    coeff_B=(-1/lambda_1)*exp(simulated_u2./lambda_1).*exp(simulated_u3./lambda_1);
    coeff_C=(1/lambda_1)*(ones(N,M)+exp(simulated_u2./lambda_1)).*exp(simulated_u3./lambda_1);
    if k<=3
        partial_g2(:,k)=(sum((coeff_A./denominator),2)/M).*stacked_X(:,k);
        partial_g3(:,k)=(sum((coeff_B./denominator),2)/M).*stacked_X(:,k);
    else
        partial_g2(:,k)=(sum((coeff_B./denominator),2)/M).*stacked_X(:,k);
        partial_g3(:,k)=(sum((coeff_C./denominator),2)/M).*stacked_X(:,k);
    end
end
partial_binding_ftn=multiplied_BZ*[partial_g2;partial_g3];
gradient=-partial_binding_ftn'*BZ'*(stacked_Y-BZ*theta_bar);
hessian=(partial_binding_ftn'*BZ')*(BZ*partial_binding_ftn);
gradient_hessian=[gradient,hessian];
end

