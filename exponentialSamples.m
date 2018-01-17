function [S,gammas]=exponentialSamples(A,n,varargin)
% exponentialSamples Compute multiresolution ensemble using exponential sampling.
%
% Syntax
%__________________________________________________________________________
%
%   [S,gammas]=exponentialSamples(A,n)
%
%   [S,gammas]=exponentialSamples(__,Name,Value)
%
%
% Description
%__________________________________________________________________________
%
%   [S,gammas]=exponentialSamples(A,n) computes exponential sampling 
%       ensemble with 'n' partitions.
%
%   [S,gammas]=exponentialSamples(__,Name,Value) additionally customizes the
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
%   'GammaMinSamples' -- Number of samples to use for 'gamma_min' search.
%                      10 (default)| scalar
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
% See Also hierarchicalConsensus, eventSamples

% Version: 1.1-alpha1
% Date: Tue Jan 16 18:15:02 EST 2018
% Author: Lucas Jeub
% Email: ljeub@iu.edu

parseArgs=inputParser();
checkFunction=@(x) isa(x,'function_handle');
addParameter(parseArgs,'Optimizer',@(B) iterated_genlouvain(B,[],0,1,...
    'moverandw'),checkFunction);
addParameter(parseArgs,'Modularity',@modularity,checkFunction)
addParameter(parseArgs,'GammaMinSamples',10,@(x) isnumeric(x) && isscalar(x))

parse(parseArgs,varargin{:});


mod_fun=parseArgs.Results.Modularity;
optimizer=parseArgs.Results.Optimizer;
samples=parseArgs.Results.GammaMinSamples;

[gamma_min, gamma_max]=gammaRange(A,'modularity',mod_fun,...
    'Samples',samples);

gammas=logspace(log10(gamma_min),log10(gamma_max),n);
% compute partitions
parfor i=1:n
    S(:,i)=optimizer(mod_fun(A,gammas(i)));
end

end


