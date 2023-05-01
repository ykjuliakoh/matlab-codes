function value=CV(h)
global X n

first=0;
second=0;

for i =1:n
    X_i=X(i,1);
    for j=1:n
        X_j=X(j,1);
        first=first+gauss_ker_con(abs((X_j-X_i)/h));
        second=second+gauss_ker((X_i-X_j)/h);
    end
end
value=(1/(n^2*h)*first)-2/(n*(n-1)*h)*(second-n*gauss_ker(0));
