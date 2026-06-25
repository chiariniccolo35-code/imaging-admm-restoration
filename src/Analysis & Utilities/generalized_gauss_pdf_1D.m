function pd = generalized_gauss_pdf_1D(m,s,q,x)
% input:
%      m:  mean
%      s:  standard deviation
%      q:  shape parameter
%      x:  x value
% output:
%   pd: p(x)    

% compute the GG pdf shape parameter
alpha = s * sqrt( gamma(1/q) / gamma(3/q) );
% compute the probability density values
pd    = q * exp( -(abs(x - m) / alpha).^q ) / ( 2 * alpha * gamma(1/q) ); 