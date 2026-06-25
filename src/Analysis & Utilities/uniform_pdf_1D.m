function pd = uniform_pdf_1D(m,s,x)
% input:
%   m:  mean
%   s:  standard deviation
%   x:  x value
% output:
%   pd: p(x) 

% compute the limits of the support interval
supp_x_min = m - sqrt(3) * s;
supp_x_max = m + sqrt(3) * s;

% compute the probability density values
pd = zeros(size(x));
pd( (x >= supp_x_min) & (x <= supp_x_max) ) = 1 / (supp_x_max - supp_x_min);