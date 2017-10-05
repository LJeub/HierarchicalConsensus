function [Sall,p]=allPartitions(Sc,Tree)
% allPartitions Find all cuts of a consensus tree
%
% Syntax
%__________________________________________________________________________
%
%   [Sall,p]=allPartitions(Sc,Tree)
%
% Description
%__________________________________________________________________________
%
%   [Sall,p]=allPartitions(Sc,Tree) computes all cuts of a hierarchical
%       tree (e.g. based the output of hierarchicalConsensus). 
%
% Input Arguments
%__________________________________________________________________________
%
%   Sc -- Finest level partition for the tree given as a vector of cluster
%         assignments
%
%   Tree -- Tree indicating how to merge clusters in 'Sc' to construct
%           larger-scale clusters. 'Tree' is an edge list, where
%           'Tree(i,1)' is the index of the coarser cluster, 'Tree(i,2)' is
%           the index of the finer cluster to be merged and 'Tree(i,3)' is
%           a measure of similarity between the clusters to be merged. (In
%           the case of a tree obtained from 'hierarchicalConsensus' this
%           is the average value of the null matrix used to split the
%           clusters.
%
% Output Arguments
%__________________________________________________________________________
%
%   Sall -- All partitions obtained by cuts of the hierarchy given as a
%           matrix
%
%   p -- The values of the similarity between clusters merged in Sall
%
% See also hierarchicalConsensus

% Version: 1.0
% Date: Thu Oct  5 17:09:19 EDT 2017
% Author: Lucas Jeub
% Email: ljeub@iu.edu

if isempty(Tree)
    Sall=Sc;
    p=0;
else
[~,s]=sort(Tree(:,3),'descend');
Tree=Tree(s,:);
Sall=Sc;
p(1)=Tree(1,3);
for i=1:size(Tree,1)
    if ~isequal(Tree(i,3),p(end))
        Sall(:,end+1)=Sc;
        p(end+1)=Tree(i,3);
    end
    Sc(Sc==Tree(i,2))=Tree(i,1);
end
end
Sall(:,end+1)=Sc;
p(end+1)=0;
end
    
