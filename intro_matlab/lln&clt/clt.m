rng(2);
%population mean=0.5, population variance=0.25
%N=10
x1=zeros(10,1000);
for i=1:1000
    x1(:,i)=rand(10,1);
end
x_1=(x1>=0.5);
mean_x1=mean(x_1,1);
stan_x1=zeros(1,1000);
for i=1:1000
    stan_x1(1,i)=sqrt(10)*(mean_x1(1,i)-0.5)/sqrt(0.25);
end

[n1,x1]=hist(stan_x1, 50);
bar(x1, n1/sum(n1));
xlabel('Bin Locations')
ylabel('Relative Frequency')

N=100;
x2=[];
for i=1:1000
    x2(:,i)=rand(100,1);
end
x_2=(x2>=0.5);
mean_x2=mean(x_2,1);
stan_x2=zeros(1,1000);
for i=1:1000
    stan_x2(1,i)=sqrt(100)*(mean_x2(1,i)-0.5)/sqrt(0.25);
end


[n2,x2]=hist(stan_x2, 50);
bar(x2, n2/sum(n2));
xlabel('Bin Locations')
ylabel('Relative Frequency')

N=1000;
x3=zeros(1000,1000);
for i=1:1000
    x3(:,i)=rand(1000,1);
end
x_3=(x3>=0.5);
mean_x3=mean(x_3,1);
stan_x3=zeros(1,1000);
for i=1:1000
    stan_x3(1,i)=sqrt(1000)*(mean_x3(1,i)-0.5)/sqrt(0.25);
end


[n3,x3]=hist(stan_x3, 50);
bar(x3, n3/sum(n3));
xlabel('Bin Locations')
ylabel('Relative Frequency')
