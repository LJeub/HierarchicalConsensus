function [Sall,thresholds]=allPartitions(Sc,Tree)
% allPartitions Find all cuts of a consensus tree
%
% Syntax
%__________________________________________________________________________
%
%   [Sall,thresholds]=allPartitions(Sc,Tree)
%
% Description
%__________________________________________________________________________
%
%   [Sall,thresholds]=allPartitions(Sc,Tree) computes all cuts of a hierarchical
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
%           a measure of similarity between the clusters to be merged.
%
% Output Arguments
%__________________________________________________________________________
%
%   Sall -- All partitions obtained by cuts of the hierarchy given as a
%           matrix
%
%   thresholds -- The values of the similarity between clusters merged in Sall
%
% See also hierarchicalConsensus

% Version: 1.1-alpha1
% Date: Tue Jan 16 18:15:01 EST 2018
% Author: Lucas Jeub
% Email: ljeub@iu.edu

if isempty(Tree)
    Sall=Sc;
    thresholds=0;
else
    Tree=sort_tree(Tree);
    Sall=Sc;
    thresholds(1)=Tree(1,3);
    com=Tree(1,1);
    for i=1:size(Tree,1)-1
        Sc=merge(Sc,Tree(i:end,:));
        skipped = isequal(Sc,Sall(:,end));
        if ~isequal(Tree(i+1,1),com) && ~skipped
            if thresholds(end)>Tree(i+1,3)
                Sall(:,end+1)=Sc;
                thresholds(end+1)=Tree(i+1,3);
            else
                warning('Similarity weakly inconsistent with hierarchy, a level has been collapsed')
            end
            com=Tree(i+1,1);
        end
        if skipped
            thresholds(end)=Tree(i+1,3);
        end
    end
    Sc=merge(Sc,Tree(end,:));
    Sall(:,end+1)=Sc;
    thresholds(end+1)=0;
end

end

function Sc=merge(Sc, Tree)
ind=find(Sc==Tree(1,2));
if isempty(ind)
    merge_ahead=find(Tree(:,1)==Tree(1,2));
    if ~isempty(merge_ahead)
        warning('Similarity inconsistent with hierarchy, a level has been collapsed');
        for j=1:length(merge_ahead)
            Sc=merge(Sc,Tree(merge_ahead(j):end,:));
        end
        ind=find(Sc==Tree(1,2));
    end
end
Sc(ind)=Tree(1,1);
end

function Tree=sort_tree(Tree)
% make sure that merges to the same group have the same distance value
[g,~,e]=unique(Tree(:,1));
G=sparse(1:size(Tree,1),e,true);
for i=1:size(G,2)
    if length(unique(Tree(G(:,i),3)))>1
        warning('merge distance not consistent for group %u and has been fixed by averaging',g(i))
        Tree(G(:,i),3)=mean(Tree(G(:,i),3));
    end
end
% sort tree (make sure merges to the same group are consecutive)
[~,s]=sortrows(Tree(:,[3,1]));
Tree=Tree(s(end:-1:1),:);
end


