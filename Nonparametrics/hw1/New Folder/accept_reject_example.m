
% F(x)=3+x^2/3-x^3/3
% f(x)=1+2x/3-x^2
% max=10/9
global X n x_made



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

hist(x_made,30)
title('histogram of generated x');

sig_hat=std(x_made);
h_r=1.06*sig_hat*N^(-0.2); % nu=2 %


% example x0=0.3 %
X=x_made;
x0=0:0.01:1;
f_hat_x0=zeros(length(x0),1);

for j=1:length(x0)
    x0_j=x0(j);
    for i=1:N
      f_hat_x0(i+1,1)=f_hat_x0(i,1)+gauss_ker((X(i,1)-x0_j)/h_r)/(N*h_r);
    end
end

x=x0';
f_x0_true=accept_reject_f(x0);

n=length(x);

