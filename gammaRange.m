function [gamma_min,gamma_max]=gammaRange(A,accuracy,varargin)
% gammaRange Compute range of gamma values.
%
% Syntax
%__________________________________________________________________________
%
%   [gamma_min,gamma_max]=gammaRange(A,accuracy)
%
%   [gamma_min,gamma_max]=gammaRange(__,Name,Value)
%
%
% Description
%__________________________________________________________________________
%
%   [gamma_min,gamma_max]=gammaRange(A,accuracy) computes range of 'gamma'
%       values with given accuracy.
%
%   [gamma_min,gamma_max]=gammaRange(__,Name,Value) additionally customizes 
%       the behavior of the function by e.g. using a different algorithm to
%       optimize modularity or using a different modularity-like quality
%       function.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   A -- Adjacency matrix of the network
%
%   accuracy -- Accuracy for the bisection search to compute 'gamma_min'.
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
%   'GammaMinBound' -- Lower bound for 'gamma_min' search. 
%                      0 (default)| scalar
%
%
% Output Arguments
%__________________________________________________________________________
%
%   gamma_min -- Smallest value of 'gamma' for which the network splits
%                into communities.
%
%   gamma_max -- Smallest value of 'gamma' for which the network is split
%                into singleton communities.
%
% See Also eventSamples, exponentialSamples, hierarchicalConsensus

% Version: 1.0
% Date: Thu Oct  5 17:09:19 EDT 2017
% Author: Lucas Jeub
% Email: ljeub@iu.edu

parseArgs=inputParser();
checkFunction=@(x) isa(x,'function_handle');
addParameter(parseArgs,'Optimizer',@(B) iterated_genlouvain(B,[],0,1,...
    'moverandw'),checkFunction);
addParameter(parseArgs,'Modularity',@modularity,checkFunction)
addParameter(parseArgs,'GammaMinBound',0,@(x) isnumeric(x) && isscalar(x))

parse(parseArgs,varargin{:});

mod_fun=parseArgs.Results.Modularity;
optimizer=parseArgs.Results.Optimizer;
bound=parseArgs.Results.GammaMinBound;

A=sparse(A);
nc_min=nComponents(A);
n=length(A);
% find smallest value of gamma necessary to make all interactions negative
gamma_max=max(max(div_0((A+A')/2,((A+A')/2-modularity(A,1)))));

% bisection search to find gamma_min
a=bound;
b=1;
check=@(gamma) max(optimizer(mod_fun(A,gamma)))>nc_min;
while ~check(b)
    b=2*b;
end
if check(a)
    gamma_min=a;
else
    while b-a>accuracy
        mid=(a+b)/2;
        if check(mid)
            b=mid;
        else
            a=mid;
        end
    end
    gamma_min=b;
end

% avoid sparse output
gamma_min=full(gamma_min);
gamma_max=full(gamma_max);
end

function c=nComponents(A)
% use dmperm to find the number of weakly connected components
[~,~,r,~]=dmperm(spones(A+A')+speye(length(A)));
c=length(r)-1;
end

function A=div_0(A,B)
% DIV_0 pointwise division such that 0/0=0
ind=find(A);
A(ind)=A(ind)./B(ind);
end
