global y_r x_r n

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1_new=x1_big(1:500,:);
x2_new=x2_big(1:500,:);
y_new=y_big(1:500,:);
v_new=v_big(1:500,:);

theta_hat_big=[];
    %[ theta1=theta(1), theta2=theta(2)]%
initial_theta_row=[1; -0.5]';

for i=1:1000
    x1_r=x1_new(:,i);
    x2_r=x2_new(:,i);
    x_r=[x1_r x2_r];
    [n,k]=size(x_r);
    y_r=y_new(:,i);
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
save result_n_500_probit.dat theta_hat_big -ascii;


