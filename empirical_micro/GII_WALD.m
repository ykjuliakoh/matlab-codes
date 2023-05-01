
function gradient_hessian=GII_WALD(X,Y,Z,theta_hat,epsilon,beta)
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

A=zeros(p,p*(M+1));
for m=1:M+1
    if m==1
        A(:,1:p*m)=eye(p);
    else
        A(:,p*(m-1)+1:p*m)=-eye(p)/M;
    end
end
S=zeros(N,p*(M+1));
for i=1:N
    Z_i=Z(i,:);
    smoothed_u2_i=smoothed_u2(i,:);smoothed_u3_i=smoothed_u3(i,:);
    for m=1:M+1
        if m==1
            S(i,1:p)=(stacked_Y(i,1)-BZ(i,:)*theta_hat)*BZ(i,:);
        else
            S(i,p*(m-1)+1:p*m)=[(smoothed_u2_i(1,m-1)-Z_i*theta_simul(1:11,m-1))*Z_i,(smoothed_u3_i(1,m-1)-Z_i*theta_simul(12:22,m-1))*Z_i];
        end
    end
end
H=BZ'*BZ;
initial_S_sum=zeros(p*(M+1),p*(M+1));
for i=1:N
    S_i=S(i,:);
    S_sum=initial_S_sum+S_i'*S_i;
    initial_S_sum=S_sum;
end
V=A*(S_sum./N)*A';
W=inv(inv(H)*V*inv(H));

gradient=partial_binding_ftn'*(theta_bar-theta_hat);
hessian=partial_binding_ftn'*partial_binding_ftn;
gradient_hessian=[gradient,hessian];
end


