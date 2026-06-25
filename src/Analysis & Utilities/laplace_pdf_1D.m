function pd = laplace_pdf_1D(m,s,x)
% input:
%   m:  mean
%   s:  standard deviation
%   x:  x value
% output:
%   pd: p(x)    

% compute the probability density values
pd = exp( -sqrt(2) * abs(x - m) / s ) / (sqrt(2) * s);