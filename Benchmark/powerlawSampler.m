function x=powerlawSampler(n,t,x_min,x_max)
% powerlawSampler Sample values from a truncated powerlaw distribution
%
% Syntax
%__________________________________________________________________________
%
%   x=powerlawSampler(n,t,x_min,x_max)
%
%
% Description
%__________________________________________________________________________
%
%   x=powerlawSampler(n,t,x_min,x_max) samples 'n' values from a truncated
%       powerlaw distribution with exponent 't', minimum cutoff 'x_min' and
%       maximum cutoff 'x_max'.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   n -- number of samples
%
%   t -- exponent
%
%   x_min -- minimum cutoff
%
%   x_max -- maximum cutoff
%
%
% Output Arguments
%__________________________________________________________________________
%
%   x -- sampled values
%
%
% See also hierachicalBenchmark

% Version: 1.0
% Date: Thu Oct  5 17:09:19 EDT 2017
% Author: Lucas Jeub
% Email: ljeub@iu.edu

y=rand(n,1);
if t~=-1
    x=((x_max^(t+1)-x_min^(t+1))*y+x_min^(t+1)).^(1/(t+1));
else
    x=x_min*(x_max/x_min).^y;
end

end
