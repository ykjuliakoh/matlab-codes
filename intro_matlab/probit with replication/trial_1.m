global y_r x_r N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1_big=load('data_x1_big.dat');
x2_big=load('data_x2_big.dat');
v_big=load('data_v_big.dat');

y_star_big=1*x1_big+(-0.5)*x2_big+v_big;

y_big=(y_star_big >= 0);

theta_hat_big=[];
initial_theta_row=[1; -0.5]';
N=length(x1_big);
for i=1:1000
    x1_r=x1_big(:,i);
    x2_r=x2_big(:,i);
    x_r=[x1_r x2_r];
    y_r=y_big(:,i);
      format long; warning off all; optimset('Tolfun',1e-20,'Tolx',1e-18,'MaxFunEvals',300*length(initial_theta_row));
      [theta_1, obj_1,flag_1]=fminsearch(@sample_obj_probit,initial_theta_row); 
      if flag_1==1
        result=[theta_1, obj_1,flag_1];  
      else
        theta_0=theta_1;
        format long; warning off all; optimset('Tolfun',1e-20,'Tolx',1e-18,'MaxFunEvals',300*length(initial_theta_row));
        [theta_2, obj_2,flag_2]=fminsearch(@sample_obj_probit,theta0); 
        result=[theta_2, obj_2,flag_2];
      end
          theta_hat_simplex=result;
          theta_hat_big(i,:)=theta_hat_simplex;
end

    
