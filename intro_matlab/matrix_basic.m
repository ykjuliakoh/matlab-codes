
%% MATRIX CONSTRUCTION %%

A = eye(5);% identity matrix
B = [1 2];%(1 by 2)
B_trans = B';
C = [3;4];%(2 by 1)
D = [1 2; 3 4];
E = [5 6; 7 9];

%Matrix multiplication
BC = B * C;
DE = D * E; %(2 by 2)
DE_e = D .* E;%(2 by 2)

%determinant and eigenvalues of matrix
det_D = det(D);
eig_D = eig(D);
inv_D = inv(D);
det_E = det(E);eig_E = eig(E);inv_E = inv(E);
F = [1 1 ; 1 1];det_F = det(F); eig_F = eig(F);