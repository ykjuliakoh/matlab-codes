function value=CV(h)
global X n

for i =1:n
    for j=1:n
        first(i,j)=gauss_ker_con((X(i,1)-X(j,1))/h);
        second(i,j)=gauss_ker((X(i,1)-X(j,1))/h);
    end
end
value=1/(n^2*h)*sum(first(:))+1/(n*(n-1)*h)*(sum(second(:))-n*gauss_ker(0));
