rng(1);

x=0.8;

u_s=randn(1000,1);
u_s_x=u_s-x;

v=[];

for i=1:1000
    v(i,1)=gauss_ker(u_s_x(i,1));
end

est_g_conv=sum(v)/1000;

g_conv=gauss_ker_con(x);