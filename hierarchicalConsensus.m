function [Sc,Tree]=hierarchicalConsensus(Ssample,varargin)
% HierarchicalConsensus Compute hierarchical consensus partition.
%
% Syntax
%__________________________________________________________________________
%
%   [Sc,Tree]=hierarchicalConsensus(Ssample)
%
%   [Sc,Tree]=hierarchicalConsensus(Ssample,alpha)
%
%   [Sc,Tree]=hierarchicalConsensus(__,Name,Value)
%
%
%
% Description
%__________________________________________________________________________
%
%   [Sc,Tree]=hierarchicalConsensus(Ssample) returns the final consensus
%       partition and consensus tree for significance level 'alpha=0.05'.
%
%   [Sc,Tree]=hierarchicalConsensus(Ssample,alpha) returns the final
%       consensus partition and consensus tree for significance level
%       'alpha'.
%
%   [Sc,Tree]=hierarchicalConsensus(__,Name,Value) additionally customizes
%       the behavior of the function using 'Name','Value' pair arguments.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   Ssample -- Ensemble of input partitions given as a matrix where each
%              column encodes a partition.
%              matrix
%
%   alpha -- Significance level for testing whether a pair of nodes is
%            significantly less frequently co-assigned than would be
%            expected by chance.
%            0 <= alpha <= 1 (default: alpha=0.05)
%
%
% Name-Value Pair Arguments
%__________________________________________________________________________
%
% Parameter names can be abbreviated and are not case sensitive.
%
%   'Iterations' -- Number of iterations to use for identifying the
%                   consensus partition at each step.
%                   size(Ssample,2) (default) | integer
%
%   'CoclassificationMatrix' -- Precomputed matrix of node-pair
%                               coclassifications.
%                               optional, useful to avoid unnecessary
%                               computation.
%
%   'Verbose' -- Display output as the algorithm progresses.
%                true (default) | false
%
%   'NullModel' -- Function handle to null model to use for the
%                  co-assignment matrix.
%                  @localPermModel (default) | @permModel | ...
%
%   'PoissBinApprox' -- Function for approximating the quantiles of the
%                       Poisson-Binomial distribution.
%                       @normApprox (default) | @sampleApprox | ...
%
%   'Optimizer' -- Function for finding "optimal" partitions for a
%                  modularity-like quality function with modularity matrix
%                  'B'
%                  @(B) iterated_genlouvain(B,[],0,1,'moverandw') (default) |
%                  function handle
%
%   'SimilarityFunction' -- Specify function used to summarize entries of
%                           the coclassification matrix into a single
%                           similarity value:
%
%                               'mean' (default) -- use mean value
%
%                               'min' -- use minimum value
%
%                               function handle -- supply a custom function
%                           
%
%   'SimilarityType' -- Specify which entries of the coclassification
%                       matrix to use to compute the similarity value:
%
%                           'all' (default) -- use all entries for nodes in
%                                              the cluster
%
%                           'linkage' -- use only entries between nodes in
%                                        different subclusters, e.g. only 
%                                        the entries marked with 'y' in the
%                                        example below:
%
%                                           1 | x x x y y y y y
%                                           1 | x x x y y y y y
%                                           1 | x x x y y y y y
%                                           2 | y y y x x y y y
%                                           2 | y y y x x y y y
%                                           3 | y y y y y x x x
%                                           3 | y y y y y x x x
%                                           3 | y y y y y x x x
%
%
% Output Arguments
%__________________________________________________________________________
%
%   Sc -- Finest level partition of the consensus hierarchy
%         vector
%
%   Tree -- Hierarchical tree indicating how to merge clusters in Sc to
%           reconstruct coarser clusters. Tree is a matrix where each row
%           represents an edge of the tree. The first element of a row
%           gives the index of the coarse cluster that the finer cluster
%           given by the second element is merged into.
%           matrix (m x 2)
%
% See Also eventSamples, drawHierarchy, consensusPlot, coclassificationMatrix,
% allPartitions, localPermModel, permModel, normApprox, sampleApprox

% Version: 1.1.1
% Date: Thu  8 Mar 2018 15:34:46 CET
% Author: Lucas Jeub
% Email: ljeub@iu.edu

[N,L]=size(Ssample);

parseArgs=inputParser();
checkFunction=@(x) isa(x,'function_handle');
addRequired(parseArgs,'sample',@ismatrix);
addOptional(parseArgs,'alpha',0.05,@(x) isnumeric(x)&&(0<=x)&&(x<=1));
addParameter(parseArgs,'NullModel',@localPermModel,checkFunction);
addParameter(parseArgs,'PoissBinApprox',@normApprox,checkFunction);
addParameter(parseArgs,'Optimizer',@(B) iterated_genlouvain(B,[],0,1,...
    'moverandw'),checkFunction);
addParameter(parseArgs,'Iterations',L,@(x) isnumeric(x) && isscalar(x));
addParameter(parseArgs,'Verbose',true,@islogical);
addParameter(parseArgs,'CoclassificationMatrix',[],@(x) ismatrix(x))
addParameter(parseArgs,'SimilarityFunction','mean');
addParameter(parseArgs,'SimilarityType','all',@(x) ismember(x,{'all','linkage'}))

parse(parseArgs,Ssample,varargin{:});

NUM_TOL=10^-10;
% constant used in genlouvain for selecting positive moves. This value is
% added to the null model to avoid splitting a pair of nodes that is
% co-classified exactly as many times as expected

% switch on or off verbose output
    function nodisp(~)
    end

if parseArgs.Results.Verbose
    mydisp=@disp;
else
    mydisp=@nodisp;
end

% build co-classification matrix unless provided
A=parseArgs.Results.CoclassificationMatrix;
if isempty(A)
    A=coclassificationMatrix(Ssample);
end

% extract remaining arguments
alpha=parseArgs.Results.alpha;
null_dist=parseArgs.Results.NullModel;
poiss_bin_approx=parseArgs.Results.PoissBinApprox;
n_iter=parseArgs.Results.Iterations;
optimizer=parseArgs.Results.Optimizer;
funInput=parseArgs.Results.SimilarityFunction;
if ischar(funInput)
    switch funInput
        case 'mean'
            simfun=@mean;
        case {'min','minimum'}
            simfun=@min;
        otherwise
            error('unknown Function %s',funInput);
    end
elseif isa(funInput,'function_handle')
    simfun=funInput;
else
    error('Function input must be string or function handle');
end
simtype=parseArgs.Results.SimilarityType;
delete(parseArgs);

Tree=zeros(0,3);
Sc=ones(N,1);
coms_old=0;
coms_new=0;
coms=1;
level=1;
% loop unitl no new communities are found
while coms_new~=coms
    G=sparse(1:N,Sc,true);
    coms_new=coms;
    mydisp(sprintf('level %u:',level));
    mydisp(sprintf('splitting %u new communities:',coms_new-coms_old));
    for i=coms_old+1:size(G,2)
        ind=find(G(:,i));
        if length(ind)>1
            A_it=A(ind,ind);
            p=poiss_bin_approx(alpha,null_dist(Ssample(ind,:)));
            if any(p(:)==0)
                mydisp(['Null model has 0-valued entries,'...
                    'increase ensemble size or reduce significance level '...
                    'for meaningful results']);
            end
            B=A_it-p+2*NUM_TOL;
            while true
                if any(B(:)<=0)
                    Sc_it=zeros(length(B),n_iter);
                    parfor c=1:n_iter
                        Sc_it(:,c)=optimizer(B);
                    end
                    A_it=coclassificationMatrix(Sc_it);
                    if isequal(A_it,double(logical(A_it)))
                        Sc_it=Sc_it(:,1);
                        break;
                    else
                        B=A_it-poiss_bin_approx(alpha,null_dist(Sc_it))+2*NUM_TOL;
                    end
                else
                    Sc_it=ones(length(B),1);
                    break;
                end
            end
            
            if max(Sc_it)>1
                Sc(ind)=Sc_it+coms;
                c_it=max(Sc_it);
                Cit=A(ind,ind);
                switch simtype
                    case 'linkage'
                        mask=bsxfun(@ne,Sc_it,Sc_it');
                        val=Cit(mask);
                    case 'all'
                        val=Cit(:);
                    otherwise
                        error('Unknown similarity type');
                end
                val=simfun(val);
                if any(val<=Tree(Tree(:,2)==i,3))
                    warning('Similarity function is not consistent with Tree for group %u',i)
                end

                % update tree
                Tree(end+(1:c_it),:)=[repmat(i,c_it,1),coms+(1:c_it)',repmat(val,c_it,1)];
                coms=coms+c_it;
                mydisp(sprintf('split com %u into %u coms',i-coms_old,max(Sc_it)))
            end
        end
    end
    level=level+1;
    coms_old=coms_new;
end

end


