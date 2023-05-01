
rng(1);

m=1000; %Replication number%

x1_big=normrnd(0,1,[1000,m]);
x2_big=exprnd(1/2, [1000,m]);
v_big=normrnd(0,1,[1000,m]);

save data_x1_big.dat x1_big -ascii;
save data_x2_big.dat x2_big -ascii;
save data_v_big.dat v_big -ascii;