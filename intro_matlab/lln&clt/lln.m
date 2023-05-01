rng(1);

x1=rand(10,1);
x_1=(x1<=0.5);
mean_10=mean(x_1);
%0.8

x2=rand(100,1);
x_2=(x2<=0.5);
mean_100=mean(x_2);
%0.46

x3=rand(1000,1);
x_3=(x3<=0.5);
mean_1000=mean(x_3);
%0.489

x4=rand(10000,1);
x_4=(x4<=0.5);
mean_10000=mean(x_4);
%0.4981

mean_lln=[mean_10; mean_100; mean_1000; mean_10000];

