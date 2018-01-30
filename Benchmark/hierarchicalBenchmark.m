function [A,S]=hierarchicalBenchmark(n,p,varargin)
% hierarchicalBenchmark Hierarchical community strucuture benchmark
%
% Syntax
%__________________________________________________________________________
%
%   [A,S]=hierarchicalBenchmark(n,p)
%
%   [A,S]=hierarchicalBenchmark(__,Name,Value)
%
%
% Description
%__________________________________________________________________________
%
%   [A,S]=hierarchicalBenchmark(n,p) samples adjacency matrix and
%       hierarchical community structure based on a degree-corrected
%       stochastic block model with 'n' nodes.
%
%   [A,S]=hierarchicalBenchmark(__,Name,Value) allows to customize
%       parameters like the degree distribution, number and size
%       distribution of communities.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   n -- number of nodes
%
%   p -- vector of probabilities (sum(p)==1) specifying the fraction of
%        edges allocated to each level of the hierarchical structure. The
%        number of hierarchical levels generated is length(p)-1 and a
%        fraction p(1) of edges is not constraint by any community
%        structure.
%
%
% Name-Value Pair Arguments
%__________________________________________________________________________
%
%   'SplitDistribution' -- Function handle or cell array of size 
%                          'length(p)-1' of function handles, that returns 
%                          (possibly random) number of communities a 
%                          community should be split into at the next level
%                          of the hierarchy. Use a cell array to specify
%                          different distributions for each hierarchical
%                          level.
%                          @() max(2,poissrnd(4)) (default) | function
%                          handle | cell array
%
%   'SizeDistribution' -- Function handle or cell array of size 
%                         'length(p)-1' of function handles. The function(s) 
%                         takes the number of communities as input and 
%                         returns a vector of probabilities for a node to
%                         be assigned to a given community. Use a cell
%                         array to specify different distributions for each
%                         hierarchical level.
%                         @(nc) dirichletSampler(1.5,nc) (default) | 
%                         function handle | cell array
%
%   'DegreeDistribution' --  Function handle specifying the distribution of
%                            expected degrees for the network. The function
%                            takes the number of nodes as input and returns
%                            a sequence of expected degrees.
%                            @(n) powerlawSampler(n,-2,5,70) (default) |
%                            function handle | cell array
%
%
% Output Arguments
%__________________________________________________________________________
%
%   A -- Adjacency matrix
%
%   S -- Planted partitions as a matrix where each column corresponds to a
%        hierarchical level. A fraction p(i) edges is constraint to lie 
%        within groups in partition S(:,i). S(:,1) is always the all-ones 
%        vector.
%
%
% See also dirchletSampler, powerlawSampler, blockModelSampler

% Version: 1.1
% Date: Tue 30 Jan 2018 18:22:37 EST
% Author: Lucas Jeub
% Email: ljeub@iu.edu

parseArgs=inputParser();
addParameter(parseArgs,'SplitDistribution',@() max(2,poissrnd(4)));
addParameter(parseArgs,'SizeDistribution',@(nc) dirichletSampler(1.5,nc));
addParameter(parseArgs,'DegreeDistribution',@(n) powerlawSampler(n,-2,5,70));
addParameter(parseArgs,'ForceConnected',true);

parse(parseArgs,varargin{:});
options=parseArgs.Results;

S=ones(n,length(p));
k=options.DegreeDistribution(n);
W=sum(k)*p(1);
A=blockModelSampler(W,S(:,1),k);
for i=1:length(p)-1
    alpha=p(i+1);
    nc=max(S(:,i));
    G=sparse(1:n,S(:,i),true,n,nc);
    coms=0;
    for c=1:nc
        sc=options.SplitDistribution();
        sc=max(round(sc),1);
        sizes=options.SizeDistribution(sc);
        S(G(:,c),i+1)=sort(discreteSample(sum(G(:,c)),sizes))+coms;
        coms=coms+sc;
    end
    A=max(A,blockModelSampler(blockMatrix(S(:,i+1),k)*alpha,S(:,i+1),k));
end

if options.ForceConnected
    A=addConnectingEdges(A,k);
end
end

function sample=discreteSample(n, weights)
weights=cumsum(weights);
weights=weights./weights(end);
sample=arrayfun(@(d) find(weights>d,1),rand(n,1));
end

function W=blockMatrix(S,k)
    G=sparse(1:length(S),S,true,length(S),max(S));
    W=diag(arrayfun(@(i) sum(k(G(:,i))),1:size(G,2)));
    W=sparse(W);
end

function A=addConnectingEdges(A,k)
    C=components(A);
    
    while(max(C)>1)
        n1=discreteSample(1,k);
        n2=find(C~=C(n1));
        n2=n2(discreteSample(1,k(n2)));
        C(C==C(n1)|C==C(n2))=min(C(n1),C(n2));
        A(n1,n2)=1;
        A(n2,n1)=1;
    end
end

function c=components(A)
    [p,~,r,~]=dmperm(spones(A+A')+speye(length(A)));
    c=zeros(length(A),1);
    for i=1:(length(r)-1)
        c(p(r(i):r(i+1)-1))=i;
    end
end
