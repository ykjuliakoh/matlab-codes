
% F(x)=3+x^2/3-x^3/3
% f(x)=1+2x/3-x^2
% max=10/9
global X n 

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

n=length(X);


%Inits = Rule of thumb value
format long; warning off all; optimset('Tolfun',1e-20,'Tolx',1e-18,'MaxFunEvals',30000*length(h_r));
[h,iter]=fminsearch(@CV,h_r); 
h_simplex=[h];

X=x_made;
x0_grid=0.1:0.01:1;
f_hat_grid=[];
f_hat_true=[];
for g=1:length(x0_grid)
    f_hat_x0=0;
   for i=1:N
       f_hat_x0=f_hat_x0+epa((X(i,1)-x0_grid(g))/h_r)/(N*h_r);
   end
   f_hat_grid(g,1)=f_hat_x0;
   f_hat_true(g,1)=accept_reject_f(x0_grid(g));
end

plot(x0_grid,f_hat_grid, x0_grid,f_hat_true')
