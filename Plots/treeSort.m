function s=treeSort(C,Sc,Tree)
% treeSort Sort coclassification matrix to be consistent with cluster tree
%
% Syntax
%__________________________________________________________________________
%
%   s=treeSort(C,Sc,Tree)
%
% Description
%__________________________________________________________________________
%
%   s=hierarchicalSort(C,Sc,Tree) sorts a coclassification matrix or other 
%       similarity matrix based on a hierarchical tree.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   C --  Coclassification matrix or other similarity matrix.
%
%   Sc -- Finest level partition for tree
%
%   Tree -- Hierarchical tree merging clusters in 'Sc'
%
%
% Output Arguments
%__________________________________________________________________________
%
%   s -- Permutation of the nodes such that C(s,s) has large entries 
%        concentrated near the diagonal
%
% See also hierarchicalSort, drawHierarchy, consensusPlot,
% hierarchicalConsensus

% Version: 1.1.1
% Date: Thu  8 Mar 2018 15:34:46 CET
% Author: Lucas Jeub
% Email: ljeub@iu.edu

% encode parents
if isempty(Tree)
    s=1:size(Sc,1);
else
A=sparse(Tree(:,2),Tree(:,1),true);
s=cell(size(A,1),1);
c=unique(Sc);
for i=1:length(c)
    ind=find(Sc==c(i));
    si=hierarchicalSort(C(ind,ind));
    s{c(i)}=ind(si);
end
for i=size(A,2):-1:1
    ind=find(A(:,i));
    if ~isempty(ind)
        Ci=zeros(length(ind),length(ind));
        for j=1:length(ind)
            for k=1:length(ind)
                Ci(j,k)=sum(sum(C(s{ind(j)},s{ind(k)})))/(length(s{ind(j)})*length(s{ind(k)}));
            end
        end
        si=hierarchicalSort(Ci);
        s{i}=vertcat(s{ind(si)});
    end
end
s=s{1};
end
end
