function obj=obj_CV_gauss(h)

global X N

sum_first=0;
sum_second=0;
for i=1:N
    for j=1:N
        sum_first=sum_first+convol_gauss((X(i,1)-X(j,1))/h);
        sum_second=sum_second+ker_gauss((X(i,1)-X(j,1))/h)*(1-(j==i));
    end
end

obj=sum_first/(N^2*h) -2* sum_second/(N*(N-1)*h);




function value=convol_gauss(t)

value=exp(-1*t^2/4)/sqrt(4*pi);


function value=ker_gauss(t)

value=exp(-0.5*t^2)/sqrt(2*pi);


