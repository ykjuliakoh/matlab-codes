function obj=obj_CV_gauss_matrix(h)

global X N

X_i_minus_X_j=zeros(N,N);

for i=1:N
        X_i_minus_X_j(:,i)=ones(N,1)*X(i,1)-X;
end

first =  sum(sum(convol_gauss(X_i_minus_X_j)))/(N^2*h);
second=  ( sum(sum(ker_gauss(X_i_minus_X_j))) -N*ker_gauss(0))/(N*(N-1)*h);

obj=first-2*second;




