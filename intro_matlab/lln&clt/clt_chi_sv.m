%%CLT using Sample Variance%%%%
%%%%%With 10 Observations%%%%%
%%%%%N=10%%%%%
variance_x_chi=var(x_chi,0,1);
sv_x_chi=[];
for i=1:1000
    sample_x_i=sample_x(:,i);
    sv_x_chi(1,i)=sqrt(10)*(mean_x(1,i)-N)/sqrt(sample_x_i);
end

[n1,xout1]=hist(sv_x_chi, 50);
bar(xout1, n1/sum(n1));
xlabel('Bin Locations')
ylabel('Relative Frequency')

%%%%%With 100 observations%%%%%%%
%%%%%%N=100%%%%%%%%
sv_x2_chi=zeros(1,1000);
for i=1:1000
    sample_x1_i=sample_x1(:,i);
    sv_x2_chi(1,i)=sqrt(100)*(mean_x1(1,i)-N)/sqrt(sample_x1_i);
end

[n2,xout2]=hist(sv_x2_chi, 50);
bar(xout2, n2/sum(n2));
xlabel('Bin Locations')
ylabel('Relative Frequency')

%%%%%clt with chi-squared with 1000 observations, using sample variance%%%%

sv_x3_chi=zeros(1,1000);
for i=1:1000
    sample_x2_i=sample_x2(:,i);
    sv_x3_chi(1,i)=sqrt(1000)*(mean_x2(1,i)-N)/sqrt(sample_x2_i);
end
[n3,xout3]=hist(sv_x3_chi, 50);
bar(xout3, n3/sum(n3));
xlabel('Bin Locations')
ylabel('Relative Frequency')