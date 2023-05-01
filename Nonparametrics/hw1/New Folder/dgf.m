
y=[];
i=1;
while i<100000
    y_i=4*randn(1,1);
    if y_i>0 && y_i<1
        y(i,1)=y_i;
        continue
    end
    i=i+1;
end

n=nnz(y) %9887
y_d=nonzeros(y);
length(y_d);
save y.dat y_d -ascii;
u=rand(length(y_d),1);
f_y=[];
for i=1:length(y_d)
    y_i=y_d(i,1);
    f_y(i,1)=1+2/3*y_i-y_i*y_i;
end
g_y=[];
for i=1:length(y_d)
    y_i=y_d(i,1);
    g_y(i,1)=4*(2/sqrt(2*pi)*exp(-(y_i^2)/2));
end

x=[];
for i=1:length(y_d)
    u_i=u(i,:);
    y_i=y_d(i,:);
    f_y_i=f_y(i,:);
    g_y_i=g_y(i,:);
    if u_i*10/9 <= f_y_i/g_y_i
        x(i,1)=y_i;
    else
        x(i,1)=0;
    end
end

x_d=nonzeros(x);

hist(x_d,50)

save x.dat x_d -ascii;    

