function value=density_ker_gauss(x0,h)

global X N

X_minus_x0_over_h= (X-ones(N,1)*x0)/h;

value=sum(ker_gauss_vector(X_minus_x0_over_h))/(N*h);

