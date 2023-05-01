global n

q_grad=[0,1,1,1]';

%%% hetero %%%
v_q_hetero=q_grad'*av_nlls_hetero*q_grad;

t_q_hetero=sqrt(n)*(theta_nlls(2)+theta_nlls(3)+theta_nlls(4))/v_q_hetero;

%%% homo %%%
v_q_homo=q_grad'*av_nlls_homo*q_grad;

t_q_homo=sqrt(n)*(theta_nlls(2)+theta_nlls(3)+theta_nlls(4))/v_q_homo;