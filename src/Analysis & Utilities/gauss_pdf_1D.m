function pd = gauss_pdf_1D(m,s,x)
% input:
%   m:  mean
%   s:  standard deviation
%   x:  x value
% output:
%   pd: p(x)    

% compute the probability density values
pd = exp( -0.5 * (((x - m) / s).^2) ) / (sqrt(2*pi) * s);