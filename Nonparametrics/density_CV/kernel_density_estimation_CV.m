
global X N

data=load('data_x.dat');
X=data(1:1000,1);
N=length(X);

sig_hat=std(X);
h_rt=1.06*sig_hat*N^(-0.2); % nu=2 %
h_min=0.5*h_rt;
h_max=3*h_rt;
grid=(h_max-h_min)/99;

h_and_CV=[];
h_grid=h_min:grid:h_max;
N_grid=length(h_grid);

for j=1:N_grid
    h_j=h_grid(j);
    h_and_CV(j,1)=h_j;
 %   h_and_CV(j,2)=obj_CV_gauss(h_j); % too time consuming
    h_and_CV(j,2)=obj_CV_gauss_matrix(h_j);
end

[h_and_CV_sort, row_h_and_CV]=sort(h_and_CV(:,2));

h_CV_gauss=h_and_CV(row_h_and_CV(1),1);

x_grid=0.1:0.01:0.99;
N_u_grid=length(x_grid);
density_gauss=[];
density_0=[];
for j=1:N_u_grid
    density_gauss(j)=density_ker_gauss(x_grid(j),h_CV_gauss);
    density_0(j)=true_f(x_grid(j));
end

plot(x_grid,density_gauss,x_grid,density_0);

