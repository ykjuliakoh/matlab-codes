%%%%%%%MSE%%%%%%%%%%%
clear
theta_results_probit=load('theta_hat_results.dat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%mean%%%%%%%%%%%%
mean(theta_results_probit(:,1)) % 1.004459842980000
mean(theta_results_probit(:,2))%-0.499880176340000

%%%%%%%var%%%%%%%%%%%%%
var(theta_results_probit(:,1)) %0.003446418329695
var(theta_results_probit(:,2))  %0.005351451322246
%%%%%%%quantile%%%%%%%%%
quantile(theta_results_probit(:,1),[0.25 0.5 0.75])
%Columns 1 through 2

 %  0.964865350000000   1.003580900000000

  %Column 3

   %1.041298800000000
quantile(theta_results_probit(:,2),[0.25 0.5 0.75])
 %Columns 1 through 2
%
 % -0.553292970000000  -0.501569155000000

  %Column 3

  %-0.448296025000000
%%%%%%%MSE%%%%%%%%%%%%%%%%
sum((theta_results_probit(:,1)-1).^2)/1000 %0.003462862110771
sum((theta_results_probit(:,2)+0.5).^2)/1000 % 0.005346114228633
%%%%%%%%HISTOGRAM%%%%%%%%%%%%
edges1=[0.7:0.017:1.2];
edges2=[-1:0.033:0];

[n1,xout1]=hist(theta_results_probit(:,1), edges1);
bar(xout1, n1/sum(n1));
title('Probit: histogram of theta1')
xlabel('Bin Locations')
ylabel('Relative Frequency')

[n2,xout2]=hist(theta_results_probit(:,2), edges2);
bar(xout2, n2/sum(n2));
title('Probit: histogram of theta2')
xlabel('Bin Locations')
ylabel('Relative Frequency')