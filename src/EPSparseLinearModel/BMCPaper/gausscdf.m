function y=gausscdf(x)
%GAUSSCDF Cumulative distribution function of standard N(0,1) Gaussian
% Y=GAUSSCDF(X) returns standard normal c.d.f. at X

y=0.5*(erf(x/sqrt(2))+1);
