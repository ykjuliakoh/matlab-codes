global h X x_made h_r
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
