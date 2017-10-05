function [S,gammas]=eventSamples(A,n,varargin)
% eventSamples Compute multiresolution ensemble using event sampling.
%
% Syntax
%__________________________________________________________________________
%
%   [S,gammas]=eventSamples(A,n)
%
%   [S,gammas]=eventSamples(__,Name,Value)
%
%
% Description
%__________________________________________________________________________
%
%   [S,gammas]=eventSamples(A,n) computes event sampling ensemble with 'n'
%       partitions.
%
%   [S,gammas]=eventSamples(__,Name,Value) additionally customizes the
%       behavior of the function by e.g. using a different algorithm to
%       optimize modularity or using a different modularity-like quality
%       function.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   A -- Adjacency matrix of the network
%
%   n -- Number of partitions to be generated
%
%
% Name-Value Pair Arguments
%__________________________________________________________________________
%
% Parameter names can be abbreviated and are not case sensitive.
%
%   'Optimizer' -- Function for finding "optimal" partitions for a
%                  modularity-like quality function with modularity matrix
%                  'B'.
%                  @(B) iterated_genlouvain(B,[],0,1,'moverandw') (default) |
%                  function handle
%
%   'Modularity' -- Function to compute modularity matrix of the form
%                   'B=modularity(A,gamma)' where 'A' is an adjacency
%                   matrix and 'gamma' is a resolution parameter. Note that
%                   the function assumes that 'B=(A+A')/2-gamma*P' for some matrix
%                   'P'.
%                   @modularity (default) | function handle
%
%   'Accuracy' -- Error tolerance for the bisection search used to find
%                 'gamma_min' the smallest value of 'gamma' for which the
%                 network is first partitioned into communities.
%                 1/n (default) | scalar
%
%   'GammaMinBound' -- Lower bound for 'gamma_min' search.
%                      0 (default)| scalar
%
%
% Output Arguments
%__________________________________________________________________________
%
%   S -- Ensemble of partitions. This is a matrix where each column
%        corresponds to a partition for a given value of 'gamma'.
%
%   gammas -- Values of 'gamma' corresponding to each partition in 'S'.
%
% See Also hierarchicalConsensus, exponentialSamples

% Version: 1.0
% Date: Thu Oct  5 17:09:19 EDT 2017
% Author: Lucas Jeub
% Email: ljeub@iu.edu

parseArgs=inputParser();
checkFunction=@(x) isa(x,'function_handle');
addParameter(parseArgs,'Optimizer',@(B) iterated_genlouvain(B,[],0,1,...
    'moverandw'),checkFunction);
addParameter(parseArgs,'Modularity',@modularity,checkFunction)
addParameter(parseArgs,'Accuracy',1/n,@(x) isnumeric(x) && isscalar(x))
addParameter(parseArgs,'GammaMinBound',0,@(x) isnumeric(x) && isscalar(x))

parse(parseArgs,varargin{:});

accuracy=parseArgs.Results.Accuracy;
mod_fun=parseArgs.Results.Modularity;
optimizer=parseArgs.Results.Optimizer;
bound=parseArgs.Results.GammaMinBound;

[gamma_min, gamma_max]=gammaRange(A,accuracy,'modularity',mod_fun,...
    'GammaMinBound',bound);

P=modularity(A,1);
A=(A+A')./2;
P=A-P;
N=length(A);
for i=1:N
    A(i,i)=0;
    P(i,i)=0;
end
ind=find(A);

% get discrete events where interactions change sign
gamma_et=full(A(ind)./P(ind));
[g_sample,~,ind2]=unique(gamma_et);

PS=sum(sum(P));
AS=sum(sum(A));
Pp=zeros(length(g_sample),1);
Ap=zeros(length(g_sample),1);
for i=1:length(g_sample)
    Pp(i)=sum(P(ind(ind2>=i)));
    Ap(i)=sum(A(ind(ind2>=i)));
end

Pp=full([PS;Pp]);
Ap=full([AS;Ap]);
g_sample=[0;g_sample];
b_sample=(g_sample.*(PS-Pp)-(AS-Ap))./(g_sample.*(PS-2*Pp)+2*Ap-AS);
b_min=(gamma_min*interp1(g_sample,PS-Pp,gamma_min,'next')-interp1(g_sample,AS-Ap,gamma_min,'next'))/...
    (gamma_min*interp1(g_sample,PS-2*Pp,gamma_min,'next')+interp1(g_sample,2*Ap-AS,gamma_min,'next'));


[b_sample,b_red]=unique(b_sample);
g_sample=g_sample(b_red);
Pp=Pp(b_red);
Ap=Ap(b_red);

% avoid outputting NaN for largest of value of gamma due to numerical error
% and 'next' Interpolant not supporting extrapolation
if b_sample(end)<1
    b_sample(end+1)=1;
    Pp(end+1)=Pp(end);
    Ap(end+1)=Ap(end);
    g_sample(end+1)=g_sample(end);
end

% interpolators to handle events
Pminus=griddedInterpolant(b_sample,PS-Pp,'next');
Pplus=griddedInterpolant(b_sample,Pp,'next');
Aminus=griddedInterpolant(b_sample,AS-Ap,'next');
Aplus=griddedInterpolant(b_sample,Ap,'next');

b=linspace(b_min,1,n)';
gammas=(Aminus(b) + b.*(Aplus(b)-Aminus(b)))./...
    ((1-b).*Pminus(b)+b.*Pplus(b));

% compute partitions
parfor i=1:n
    S(:,i)=optimizer(mod_fun(A,gammas(i)));
end

end


