x=0:0.01:1;
f_x=[];
for i=1:101
    x_i=x(:,i);
    f_x(i,1)=1+x_i*2/3-x_i*x_i;
end
g_x=[];
for i=1:101
    x_i=x(:,i);
    g_x(i,1)=4*(2/sqrt(2*pi)*exp(-(x_i^2)/2));
end

plot(x',f_x,x',g_x)
legend('f_x','g_x')