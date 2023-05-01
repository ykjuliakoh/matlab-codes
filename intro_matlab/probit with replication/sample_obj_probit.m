function objective=sample_obj_probit(theta)

% theta is a row vector %
global   x_r y_r n
  
logL=0;
for j=1:n
     xtheta=x_r(j,:)*theta';
     logL=logL+ y_r(j,1)*log(normcdf(xtheta))+ (1-y_r(j,1))*log(1-normcdf(xtheta));
end

 objective=-1*logL ;