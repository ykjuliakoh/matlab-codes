
%%CLT with Sample Variance?
rng(2);
%population mean=0.5, population variance=0.25
%N=10
x1=zeros(10,1000);
for i=1:1000
    x1(:,i)=rand(10,1);
end
x_1=(x1>=0.5);
mean_x1=mean(x_1,1);
N=10;
sample_x1=[];
for i=1:1000
    mean_x1_i=mean_x1(:,i);
    x1_i=x_1(:,i);
    sample_x1(:,i)=sumsqr(x1_i-mean_x1_i)/(N-1);
end
variance_x1=var(x1,0,1);

sv_x1=zeros(1,1000);
for i=1:1000
    variance_x1_i=variance_x1(:,i);
    sv_x1(1,i)=sqrt(10)*(mean_x1(1,i)-0.5)/sqrt(variance_x1_i);
end

[n1,xout1]=hist(sv_x1, 50);
bar(xout1, n1/sum(n1));
xlabel('Bin Locations')
ylabel('Relative Frequency')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=100;
x2=[];
for i=1:1000
    x2(:,i)=rand(100,1);
end
x_2=(x2>=0.5);
mean_x2=mean(x_2,1);

variance_x2=var(x2,0,1);
sv_x2=zeros(1,1000);
for i=1:1000
    var_x2_i=variance_x2(:,i);
    sv_x2(1,i)=sqrt(100)*(mean_x2(1,i)-0.5)/sqrt(var_x2_i);
end

[n2,xout2]=hist(sv_x2, 50);
bar(xout2, n2/sum(n2));
xlabel('Bin Locations')
ylabel('Relative Frequency')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=1000;
x3=zeros(1000,1000);
for i=1:1000
    x3(:,i)=rand(1000,1);
end
x_3=(x3>=0.5);
mean_x3=mean(x_3,1);

variance_x3=var(x3,0,1);
sv_x3=zeros(1,1000);
for i=1:1000
    var_x3_i=variance_x3(:,i);
    sv_x3(1,i)=sqrt(1000)*(mean_x3(1,i)-0.5)/sqrt(var_x3_i);
end

[n3,xout3]=hist(sv_x3, 50);
bar(xout3, n3/sum(n3));
xlabel('Bin Locations')
ylabel('Relative Frequency')