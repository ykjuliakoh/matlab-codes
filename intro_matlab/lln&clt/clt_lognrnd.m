rng(1);

mu=0;
var=1;

pop_mean=exp(mu+var/2);
pop_var=exp(2*(mu+var))-exp(2*mu+var);

%%%%%%%%%%%N=10%%%%%%%%%%%%
x=zeros(10,1000);
for i=1:1000
    x(:,i)=lognrnd(0,1,[10,1]);
end

mean_x=mean(x,1);
std_x=std(x,1);

stan_x=zeros(1,1000);
for i=1:1000
    stan_x(1,i)=sqrt(10)*(mean_x(1,i)-pop_mean)/pop_var;
end

[n,x]=hist(stan_x, 50);
bar(x, n/sum(n));
xlabel('Bin Locations')
ylabel('Relative Frequency')

%%%%%N=100%%%%%%%%%%%%%%%%%%
x1=zeros(100,1000);
for i=1:1000
    x1(:,i)=lognrnd(0,1,[100,1]);
end

mean_x1=mean(x1,1);
std_x1=std(x1,1);

stan_x1=zeros(1,1000);
for i=1:1000
    stan_x1(1,i)=sqrt(100)*(mean_x1(1,i)-pop_mean)/pop_var;
end

[n1,x1]=hist(stan_x1, 50);
bar(x1, n1/sum(n1));
xlabel('Bin Locations')
ylabel('Relative Frequency')

%%%%%%N=1000%%%%%%%%%%%%%
x2=zeros(1000,1000);
for i=1:1000
    x2(:,i)=lognrnd(0,1,[1000,1]);
end

mean_x2=mean(x2,1);
std_x2=std(x2,1);

stan_x2=zeros(1,1000);
for i=1:1000
    stan_x2(1,i)=sqrt(1000)*(mean_x2(1,i)-pop_mean)/pop_var;
end

[n2,x2]=hist(stan_x2, 50);
bar(x2, n2/sum(n2));
xlabel('Bin Locations')
ylabel('Relative Frequency')