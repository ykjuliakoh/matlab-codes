

rng(1);

x1=normrnd(0,1,[1000,1]);
x2=exprnd(1/4, [1000,1]);
v=normrnd(0,1,[1000,1]);
v_tilde=normrnd(0,1, [1000,1]);

save data_x1.dat x1 -ascii;
save data_x2.dat x2 -ascii;
save data_v.dat v -ascii;