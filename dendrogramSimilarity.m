function [Tree,isConsistent]=dendrogramSimilarity(C,Sc,Tree,varargin)
% dendrogramSimilarity Compute within-cluster similarity for hierarchy
%
% This function implements different ways to compute the within-cluster
% similarity of nodes based on the coclassification matrix for each cluster
% identified by a hierarchy.
%
% Syntax
%__________________________________________________________________________
%
%   [Tree,isConsistent]=dendrogramSimilarity(C,Sc,Tree)
%
%   [Tree,isConsistent]=dendrogramSimilarity(__,Name,Value)
%
%
% Description
%__________________________________________________________________________
%
%   [Tree,isConsistent]=dendrogramSimilarity(C,Sc,Tree) computes similarity
%       for each cluster in the Tree using the mean value of the
%       coclassification matrix constrained to nodes in the cluster. Tree 
%       is a matrix where each row represents an edge of the tree. The 
%       first element of a row gives the index of the coarse cluster that 
%       the finer cluster given by the second element is merged into. The
%       function also checks whether the similarity values are consistent
%       with the hierarchy, i.e., the similarity value for a cluster should
%       be smaller than that for its sub-clusters in the hierarchy.
%
%
%   [Tree,isConsistent]=hierarchicalConsensus(__,Name,Value) allows one to 
%       customize the function used to compute similarity values.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   C -- Coclassification matrix
%   
%   Sc -- Finest-level consensus partition
%
%   Tree -- Tree indicating how to merge clusters in Sc into larger
%           clusters.
%                           
%
% Name-Value Pair Arguments
%__________________________________________________________________________
%
% Parameter names can be abbreviated and are not case sensitive.
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
%   Tree -- Hierarchical tree indicating how to merge clusters in Sc to
%           reconstruct coarser clusters. The first and second column of 
%           'Tree' are unchanged. The third column contains the new 
%           similarity values.
%
%   isConsistent -- Boolean value indicating whether the computed
%                   similarity values are consistent with the hierarchy in
%                   the sense that the similarity value of a coarser
%                   cluster should be smaller than the similarity value of
%                   any of its sub-clusters. The function also outputs a
%                   warning for any clusters that violate this condition.
%
% See Also eventSamples, drawHierarchy, consensusPlot, coclassificationMatrix,
% allPartitions, localPermModel, permModel, normApprox, sampleApprox

% Version: 1.1
% Date: Tue 30 Jan 2018 18:22:37 EST
% Author: Lucas Jeub
% Email: ljeub@iu.edu

isConsistent=true;
if isempty(Tree)
    % Don't do anything if tree is empty
    return
end

parseArgs=inputParser();
addParameter(parseArgs,'SimilarityFunction','mean');
addParameter(parseArgs,'SimilarityType','all',@(x) ismember(x,{'all','linkage'}))
parse(parseArgs,varargin{:});

funInput=parseArgs.Results.SimilarityFunction;
if ischar(funInput)
    switch funInput
        case 'mean'
            fun=@mean;
        case {'min','minimum'}
            fun=@min;
        otherwise
            error('unknown Function %s',funInput);
    end
elseif isa(funInput,'function_handle')
    fun=funInput;
else
    error('Function input must be string or function handle');
end
type=parseArgs.Results.SimilarityType;


[merges,~,index]=unique(Tree(:,1));
Tree(:,3)=0;
for i=numel(merges):-1:1
    groups_to_merge=Tree(index==i,2);
    nodes_to_merge=find(ismember(Sc,groups_to_merge));
    Cit=C(nodes_to_merge,nodes_to_merge);
    switch type
        case 'linkage'
            mask=bsxfun(@ne,Sc(nodes_to_merge),Sc(nodes_to_merge)');
            val=Cit(mask);
        case 'all'
            val=Cit(:);
        otherwise
            error('Unknown similarity type');
    end
    val=fun(val);
    if any(val>=Tree(ismember(Tree(:,1),groups_to_merge),3))
        isConsistent=false;
        warning('Similarity function is not consistent with Tree for group %u',merges(i))
    end
    Sc(nodes_to_merge)=merges(i);
    Tree(index==i,3)=val;
end

end
