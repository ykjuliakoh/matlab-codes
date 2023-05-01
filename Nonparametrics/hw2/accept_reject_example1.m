global n X
% F(x)=3+x^2/3-x^3/3
% f(x)=1+2x/3-x^2
% max=10/9

rng(1);
M=10/9;
N=5000;

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

hist(x_made,30);
title('Histogram of Generated x')

X=x_made;
n=length(X);
sig_hat=std(x_made);
h_r=1.06*sig_hat*N^(-0.2); % nu=2 %

      format long; warning off all; optimset('Tolfun',1e-20,'Tolx',1e-18,'MaxFunEvals',3000*length(h_r));
      [h_1, obj_1,flag_1]=fminsearch(@CV,h_r); 
      if flag_1==1
        result=[h_1, obj_1,flag_1];
      else
        theta_0=h_1;
        format long; warning off all; optimset('Tolfun',1e-20,'Tolx',1e-18,'MaxFunEvals',3000*length(h_r));
        [h_2, obj_2,flag_2]=fminsearch(@CV,h_r); 
        result=[h_2, obj_2,flag_2];
      end
     CV_h=result;


X=x_made;
x0_grid=0.1:0.01:1;
f_hat_grid=[];
f_hat_true=[];
for g=1:length(x0_grid)
    f_hat_x0=0;
   for i=1:N
       f_hat_x0=f_hat_x0+gauss_ker((X(i,1)-x0_grid(g))/CV_h(1,1))/(N*CV_h(1,1));
   end
   f_hat_grid(g,1)=f_hat_x0;
   f_hat_true(g,1)=accept_reject_f(x0_grid(g));
end

plot(x0_grid,f_hat_grid, x0_grid,f_hat_true')
title('plot of f_hat with CV vs. plot of f_true')
legend('f_hat with CV','f_true')