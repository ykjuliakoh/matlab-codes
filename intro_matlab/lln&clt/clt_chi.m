rng(3);

N=10;
x_chi=zeros(10,1000);
for i=1:1000
    x_chi(:,i)=chi2rnd(N, [N 1]);
end
mean_x=mean(x_chi,1);
sample_x=[];
for i=1:1000
    mean_x_i=mean_x(:,i);
    x_chi_i=x_chi(:,i);
    sample_x(:,i)=sumsqr(x_chi_i-mean_x_i)/(N-1);
end

stan_chi=zeros(1,1000);
for i=1:1000
    stan_chi(1,i)=sqrt(N)*(mean_x(1,i)-N)/sqrt(2*N);
end


[n1,x1]=hist(stan_chi, 50);
bar(x1, n1/sum(n1));
xlabel('Bin Locations')
ylabel('Relative Frequency')

N=100;
x1_chi=zeros(100,1000);
for i=1:1000
    x1_chi(:,i)=chi2rnd(N, [N 1]);
end
mean_x1=mean(x1_chi,1);
%sample variance
sample_x1=[];
for i=1:1000
    mean_x1_i=mean_x1(:,i);
    x1_chi_i=x1_chi(:,i);
    sample_x1(:,i)=sumsqr(x1_chi_i-mean_x1_i)/(N-1);
end
stan_chi1=zeros(1,1000);
for i=1:1000
    stan_chi1(1,i)=sqrt(N)*(mean_x1(1,i)-N)/sqrt(2*N);
end

[n2,x2]=hist(stan_chi1, 50);
bar(x2, n2/sum(n2));
xlabel('Bin Locations')
ylabel('Relative Frequency')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=1000;
x2_chi=zeros(1000,1000);
for i=1:1000
    x2_chi(:,i)=chi2rnd(N, [N 1]);
end
mean_x2=mean(x2_chi,1);
%sample variance%%
sample_x2=[];
for i=1:1000
    mean_x2_i=mean_x2(:,i);
    x2_chi_i=x2_chi(:,i);
    sample_x2(:,i)=sumsqr(x2_chi_i-mean_x2_i)/(N-1);
end
%%%%%%%%%%%%%%%%%%%%
stan_chi2=zeros(1,1000);
for i=1:1000
    stan_chi2(1,i)=sqrt(N)*(mean_x2(1,i)-N)/sqrt(2*N);
end

[n3,x3]=hist(stan_chi2, 50);
bar(x3, n3/sum(n3));
xlabel('Bin Locations')
ylabel('Relative Frequency')

