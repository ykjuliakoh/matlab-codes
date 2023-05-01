global x y N S v_tilde_SN

N=1000; S=1000;
rng(1);

x1=load('data_x1.dat');
x2=load('data_x2.dat');
v=load('data_v.dat');
x=[x1,x2];
v_tilde_SN=randn(S,N);
y_star=1*x1 -0.5*x2+v;
y=[];
y=(y_star>=0);


theta0=[1, -0.5];
format long; warning off all; optimset('Tolfun',1e-20,'Tolx',1e-18,'MaxFunEvals',300*length(theta0));
[theta, obj, flag]=fminsearch(@obj_SMM_example, theta0); 
if flag==1;
   result_SMM=[theta, obj, flag];
else
   theta0=theta;    
   format long; warning off all; optimset('Tolfun',1e-20,'Tolx',1e-18,'MaxFunEvals',300*length(theta0));
   [theta2, obj2, flag2]=fminsearch(@obj_SMM_example,theta0); 
   result_SMM= [theta2, obj2, flag2];
end

