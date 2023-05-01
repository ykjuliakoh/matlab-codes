clear;
x=load('x.dat');

hist(x,50)
title('histogram of x')
hold on
f_x=[];
for i=1:length(x)
    x_i=x(i,1);
    f_x(i,1)=1+2/3*x_i-x_i*x_i;
end
plot(x,f_x,'LineWidth',1.5)

%mean(x)=0.51
sv_x=[];
sv_x=(x-0.51).^2;
sum(sv_x)/(length(x)-1)
%standard deviation=0.28
%%Pick x0=0.5
%Bandwidth using silvermann constant
h=[];
h(1,1)=0.28*2.34*length(x)^(-1/5);
h(2,1)=0.28*2.78*length(x)^(-1/5);
h(3,1)=0.28*3.15*length(x)^(-1/5);
h(4,1)=0.28*1.06*length(x)^(-1/5);
h(1:4,2)=2*h(1:4,1);
h(1:4,3)=5*h(1:4,1);
x0=0.5;
f_x0=1+2*x0/3-x0*x0; %1.0833
%%%%% WITH h*1 %%%%%%%%%%%
%Epanechnikov
k1=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(1,1);
    k1=k1+epanechnikov2(u);
end
f_hat_e=k1/(length(x)*h(1,1));%1.0667

%Biweight
k2=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(2,1);
    k2=k2+biweight2(u);
end
f_hat_b=k2/(length(x)*h(2,1));%1.0705

%Triweight
k3=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(3,1);
    k3=k3+triweight2(u);
end
f_hat_t=k3/(length(x)*h(3,1)); %1.0724

%Gaussian
k4=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(4,1);
    k4=k4+gaussian2(u);
end
f_hat_g=k4/(length(x)*h(4,1)); %1.0769

%%%%% WITH h*2 %%%%%%%%%%%
%Epanechnikov
k1_2=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(1,2);
    k1_2=k1_2+epanechnikov2(u);
end
f_hat_e2=k1_2/(length(x)*h(1,2));%1.0576

%Biweight
k2_2=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(2,2);
    k2_2=k2_2+biweight2(u);
end
f_hat_b2=k2_2/(length(x)*h(2,2));%1.0566

%Triweight
k3_2=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(3,2);
    k3_2=k3_2+triweight2(u);
end
f_hat_t2=k3_2/(length(x)*h(3,2)); %1.0566

%Gaussian
k4_2=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(4,2);
    k4_2=k4_2+gaussian2(u);
end
f_hat_g2=k4_2/(length(x)*h(4,2)); %1.0565

%%%%% WITH h*5 %%%%%%%%%%%
%Epanechnikov
k1_3=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(1,3);
    k1_3=k1_3+epanechnikov2(u);
end
f_hat_e3=k1_3/(length(x)*h(1,3));%0.9366

%Biweight
k2_3=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(2,3);
    k2_3=k2_3+biweight2(u);
end
f_hat_b3=k2_3/(length(x)*h(2,3));%0.9300

%Triweight
k3_3=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(3,3);
    k3_3=k3_3+triweight2(u);
end
f_hat_t3=k3_3/(length(x)*h(3,3)); %0.9295

%Gaussian
k4_3=0;
for i=1:length(x)
    x_i=x(i,1);
    u=(x_i-x0)/h(4,3);
    k4_3=k4_3+gaussian2(u);
end
f_hat_g3=k4_3/(length(x)*h(4,3)); %0.9277

f_hat=[];
f_hat(1:4,1)=[f_hat_e;f_hat_b;f_hat_t;f_hat_g];
f_hat(1:4,2)=[f_hat_e2;f_hat_b2;f_hat_t2;f_hat_g2];
f_hat(1:4,3)=[f_hat_e3;f_hat_b3;f_hat_t3;f_hat_g3];

f_diff=[];
f_diff=f_hat-1.0833;