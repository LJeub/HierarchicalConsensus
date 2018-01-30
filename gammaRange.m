function [gamma_min,gamma_max]=gammaRange(A,varargin)
% gammaRange Compute range of gamma values. 
%
% Syntax
%__________________________________________________________________________
%
%   [gamma_min,gamma_max]=gammaRange(A)
%
%   [gamma_min,gamma_max]=gammaRange(__,Name,Value)
%
%
% Description
%__________________________________________________________________________
%
%   [gamma_min,gamma_max]=gammaRange(A) computes range of 'gamma'
%       values that result in non-trivial partitions.
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
%                   the function assumes that 'B=(A+A')/2-gamma*P' for some
%                   matrix 'P'.
%                   @modularity (default) | function handle
%
%   'InitialGuess' -- Initial value for 'gamma_min' search. 
%                      1 (default)| scalar
%
%   'Samples' -- Number of partitions to sample to test for 'gamma_min' at
%                proposed value. More samples should result in more
%                accurate values. 
%
%
% Output Arguments
%__________________________________________________________________________
%
%   gamma_min -- Upper bound for smallest value of 'gamma' for which the 
%                network splits into communities.
%
%   gamma_max -- Smallest value of 'gamma' for which the network is split
%                into singleton communities.
%
%
% Implementation
%__________________________________________________________________________
%
% 'gamma_max' is simply the largest value of gamma for which there exist
% ferromagnetic interactions in the modularity matrix and is easy to
% compute directly based on the adjacency and modularity matrix. 
%
% 'gamma_min' does not have a closed-form solution and needs to be
% approximated numerically. This function uses an iterative algorithm that
% exploits the linearity of modularity as a function of gamma for a fixed
% partition. This means that we can directly compute the minimum value of
% gamma for which a given partition is better than the trivial partition.
% The iterative algorithm proceeds by first estimating 'gamma_min' using a
% small sample of partitions. We then sample a new set of patitions using
% 'gamma=gamma_min-epsilon' (to ensure the previous partitions are strictly
% non-optimal) and use the new sample to update 'gamma_min' and repeat. The
% algorithm stops once the new sample consists only of the trivial
% partition.
%
%
% See Also eventSamples, exponentialSamples, hierarchicalConsensus

% Version: 1.1-alpha1
% Date: Tue Jan 16 18:15:02 EST 2018
% Author: Lucas Jeub
% Email: ljeub@iu.edu

parseArgs=inputParser();
checkFunction=@(x) isa(x,'function_handle');
addParameter(parseArgs,'Optimizer',@(B) iterated_genlouvain(B,[],0,1,...
    'moverandw'),checkFunction);
addParameter(parseArgs,'Modularity',@modularity,checkFunction);
addParameter(parseArgs,'InitialGuess',1,@(x) isnumeric(x) && isscalar(x));
addParameter(parseArgs,'Samples',10,@(x) isnumeric(x) && isscalar(x));

parse(parseArgs,varargin{:});

mod_fun=parseArgs.Results.Modularity;
optimizer=parseArgs.Results.Optimizer;
initial_guess=parseArgs.Results.InitialGuess;
samples=parseArgs.Results.Samples;

A=sparse(A);
% find smallest value of gamma necessary to make all interactions negative
AT=(A+A')/2;
PT=AT-mod_fun(A,1);
gamma_max=max(max(div_0(AT,PT)));
NUM_TOL=10^-9; % small constant used to ensure that already sampled partitions are strictly non-optimal at the next iteration.

% make sure initial guess is valid
B0=mod_fun(A,0);
if any(B0(:)<=0)
    S0=optimizer(B0);
    G=sparse(1:length(S0),S0,1);
    a0=trace(G'*AT*G);
    p0=trace(G'*PT*G);
else
    a0=sum(sum(AT));
    p0=sum(sum(PT));
end
gamma_min=inf;
while gamma_min>=inf
    parfor i=1:samples
        S(:,i)=optimizer(mod_fun(A,initial_guess-NUM_TOL));
    end
    gamma_min=gamma_min_bound(AT,PT,S,a0,p0);
    initial_guess=min(2*initial_guess,gamma_max-NUM_TOL);
end
% update using convex hull idea
gamma_min_new=gamma_min;
gamma_min=inf;
while gamma_min_new<gamma_min
    gamma_min=gamma_min_new;
    parfor i=1:samples
        S(:,i)=optimizer(mod_fun(A,gamma_min-NUM_TOL));
    end
    gamma_min_new=gamma_min_bound(AT,PT,S,a0,p0);
end

gamma_min=full(gamma_min);
gamma_max=full(gamma_max);
end

function A=div_0(A,B)
% DIV_0 pointwise division such that 0/0=0
ind=find(A);
A(ind)=A(ind)./B(ind);
end

function gamma_min=gamma_min_bound(A,P,S,a0,p0)
    gamma_min=inf;
    for i=1:size(S,2)
        [u,~,e]=unique(S(:,i));
        if numel(u)>1
            G=sparse(1:size(S,1),e,1,size(S,1),numel(u));
            a=trace(G'*A*G);
            p=trace(G'*P*G);
            gamma_min=min(gamma_min,(a0-a)/(p0-p));
        end
    end
end
        
