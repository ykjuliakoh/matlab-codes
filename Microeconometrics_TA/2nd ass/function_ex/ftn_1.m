%%%%EXAMPLE2%%%%%%%%%
% Compute the value of the integrand at 2*pi/3.
x = 2*pi/3;
y = myIntegrand(x)

% Compute the area under the curve from 0 to pi.
xmin = 0;
xmax = pi;
f = @myIntegrand;
a = integral(f,xmin,xmax)


%%%%%%EXAMPLE2%%%%%%%%
values = [12.7, 45.4, 98.9, 26.6, 53.1];
[ave,stdev] = stat(values)