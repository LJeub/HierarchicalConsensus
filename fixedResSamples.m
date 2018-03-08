function S=fixedResSamples(A,n,varargin)
% fixedResSamples Compute ensemble using fixed resolution.
%
% Syntax
%__________________________________________________________________________
%
%   S=fixedResSamples(A,n)
%
%   S=exponentialSamples(__,Name,Value)
%
%
% Description
%__________________________________________________________________________
%
%   S=fixedResSamples(A,n) computes ensemble with 'n' partitions at 'gamma=1'.
%
%   S=fixedResSamples(__,Name,Value) additionally customizes the
%       behavior of the function by e.g. different 'gamma', using a different
%       algorithm to optimize modularity or using a different modularity-like
%       quality function.
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
%   'Gamma' -- Resolution parameter value 
%              1 (default) | scalar
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
%
% Output Arguments
%__________________________________________________________________________
%
%   S -- Ensemble of partitions. This is a matrix where each column
%        corresponds to a partition.
%
% See Also hierarchicalConsensus, eventSamples, exponentialSamples

% Version: 1.1.1
% Date: Thu  8 Mar 2018 15:34:46 CET
% Author: Lucas Jeub
% Email: ljeub@iu.edu

parseArgs=inputParser();
checkFunction=@(x) isa(x,'function_handle');
addParameter(parseArgs,'Gamma',1,@(x) isnumeric(x) && isscalar(x));
addParameter(parseArgs,'Optimizer',@(B) iterated_genlouvain(B,[],0,1,...
    'moverandw'),checkFunction);
addParameter(parseArgs,'Modularity',@modularity,checkFunction)

parse(parseArgs,varargin{:});

gamma=parseArgs.Results.Gamma;
mod_fun=parseArgs.Results.Modularity;
optimizer=parseArgs.Results.Optimizer;

parfor i=1:n
    S(:,i)=optimizer(mod_fun(A,gamma));
end

end


