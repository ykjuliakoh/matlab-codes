
% F(x)=3+x^2/3-x^3/3
% f(x)=1+2x/3-x^2
% max=10/9

rng(1);
M=10/9;
N=1000;

u_g_draw=rand(2*N,2); % set sufficiently large N
u_draw=u_g_draw(:,1); 
x_from_g=u_g_draw(:,2); 

x_made=[];
made=0;
i=1;
while made <N
     f_i=accept_reject_f(x_from_g(i));
     if u_draw(i)<= f_i/M 
        made=made+1;
        x_made(made,1)=x_from_g(i);
     else
         made=made;
     end
    i=i+1;
end

hist(x_made, 30)
title('accept-reject method: histogram of x')

sig_hat=std(x_made);
h4=1.06*sig_hat*N^(-0.2); % nu=2 %
h1=2.34*sig_hat*N^(-0.2);%epane
h2=2.78*sig_hat*N^(-0.2);%biweight
h3=3.15*sig_hat*N^(-0.2);
% example x0=0.5 %
X=x_made;
x0=0.5;
f_hat_x0=0;
for i=1:N
    f_hat_x0=f_hat_x0+gauss_ker((X(i,1)-x0)/h4)/(N*h4);
end

f_x0_true=accept_reject_f(x0);
