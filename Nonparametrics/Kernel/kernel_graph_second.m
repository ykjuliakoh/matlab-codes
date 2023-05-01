%create sequence of numbers 
%(first:increment:second)

u=(-1:0.01:1);
n=length(u);
%Uniform Kernel
k_0_u=[];
for i=1:n
    u_i=u(1,i);
    indi=(abs(u_i)<=1);
    k_0_u(i,1)=indi/2;
end

plot(u',k_0_u);

%Epanechnikov
k_1_u=[];
for i=1:n
    u_i=u(1,i);
    indi=(abs(u_i)<=1);
    k_1_u(i,1)=3*(1-u_i^2)/4;
end

plot(k_1_u);

%Biweight
k_2_u=[];
for i=1:n
    u_i=u(1,i);
    indi=(abs(u_i)<=1);
    k_2_u(i,1)=15*((1-u_i^2)^2)/16*indi;
end

plot(k_2_u);

%Triweight
k_3_u=[];
for i=1:n
    u_i=u(1,i);
    indi=(abs(u_i)<=1);
    k_3_u(i,1)=35*((1-u_i^2)^3)/32*indi;
end

plot(k_3_u);

%Gaussian
k_inf_u=[];
for i=1:n
    u_i=u(1,i);
    k_inf_u(i,1)=exp(-(u_i^2)/2)/sqrt(2*pi)*indi;
end


figure
subplot(2,3,1)       % add first plot in 2 x 2 grid
plot(u',k_0_u)           % line plot
title('Uniform Kernel')

subplot(2,3,2)       % add second plot in 2 x 2 grid
plot(u',k_1_u)        % scatter plot
title('Epanechnikov Kernel')

subplot(2,3,3)       % add third plot in 2 x 2 grid
plot(u',k_2_u)           % stem plot
title('Biweight Kernel')

subplot(2,3,4)       % add fourth plot in 2 x 2 grid
plot(u',k_3_u)
title('Triweight Kernel')

subplot(2,3,5)
plot(u',k_inf_u)
title('Gaussian Kernel')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Graph on same plot%%%%%%%
figure
plot(u',k_0_u, u',k_1_u, u',k_2_u, u',k_3_u, u',k_inf_u)
title('Common Second-Order Kernels')
legend('Uniform', 'Epanechnikov', 'Biweight', 'Triweight', 'Gaussian')



%%%%%%%%%%%%4TH ORDER KERNELS%%%%%%%%%

%Epanechnikov%%%%%%
k_4_1=[];
for i=1:n
    u_i=u(:,i);
    k_1_i=k_1_u(i,1);
    k_4_1(i,1)=15*(1-7*(u_i^2)/3)*k_1_i/8;
end
%%Biweight%%%%%%
k_4_2=[];
for i=1:n
    u_i=u(1,i);
    k_2_i=k_2_u(i,1);
    k_4_2(i,1)=7*(1-3*u_i^2)*k_2_i/4;
end
%%%%Triweight%%%%%%
k_4_3=[];
for i=1:n
    u_i=u(1,i);
    k_3_i=k_3_u(i,1);
    k_4_3(i,1)=27*(1-11*(u_i^2)/3)*k_3_i/16;
end
%%%%%%%%Gaussian%%%%%%
k_4_inf=[];
for i=1:n
    u_i=u(1,i);
    k_inf_i=k_inf_u(i,1);
    k_4_inf(i,1)=(3-u_i^2)*k_inf_i/2;
end

figure
plot(u',k_4_1, u',k_4_2, u',k_4_3, u',k_4_inf)
title('Fourth-Order Kernels')
legend('Epanechnikov', 'Biweight', 'Triweight', 'Gaussian')

%%%%%%%%%%%%6TH ORDER KERNELS%%%%%%%%%

%Epanechnikov%%%%%%
k_6_1=[];
for i=1:n
    u_i=u(:,i);
    k_1_i=k_1_u(i,1);
    k_6_1(i,1)=175*(1-6*(u_i^2)+33*(u_i^4)/5)*k_1_i/64;
end
%%Biweight%%%%%%
k_6_2=[];
for i=1:n
    u_i=u(1,i);
    k_2_i=k_2_u(i,1);
    k_6_2(i,1)=315*(1-22*(u_i^2)/3+143*(u_i^4)/15)*k_2_i/128;
end
%%%%Triweight%%%%%%
k_6_3=[];
for i=1:n
    u_i=u(1,i);
    k_3_i=k_3_u(i,1);
    k_6_3(i,1)=297*(1-26*(u_i^2)/3+13*(u_i^4))*k_3_i/128;
end
%%%%%%%%Gaussian%%%%%%
k_6_inf=[];
for i=1:n
    u_i=u(1,i);
    k_inf_i=k_inf_u(i,1);
    k_6_inf(i,1)=(15-10*u_i^2+u_i^4)*k_inf_i/8;
end

figure
plot(u',k_6_1, u',k_6_2, u',k_6_3, u',k_6_inf)
title('Sixth-Order Kernels')
legend('Epanechnikov', 'Biweight', 'Triweight', 'Gaussian')
%%%%%%%%%%%%%

figure
subplot(2,2,1)       % add first plot in 2 x 2 grid
           % line plot
plot(u',k_0_u, u',k_1_u, u',k_2_u, u',k_3_u, u',k_inf_u)
title('Common Second-Order Kernels')


subplot(2,2,2)       % add second plot in 2 x 2 grid
plot(u',k_4_1, u',k_4_2, u',k_4_3, u',k_4_inf)
title('Fourth-Order Kernels')

subplot(2,2,3)       % add third plot in 2 x 2 grid
plot(u',k_6_1, u',k_6_2, u',k_6_3, u',k_6_inf)
title('Sixth-Order Kernels')

