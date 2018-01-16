function s=hierarchicalSort(C)
% hierarchicalSort Sort coclassification matrix based on average-linkage clustering
%
% Syntax
%__________________________________________________________________________
%
%   s=hierarchicalSort(C)
%
% Description
%__________________________________________________________________________
%
%   s=hierarchicalSort(C) sorts a coclassification matrix or other 
%       similarity matrix based on finding the optimal leaf order for the 
%       average-linkage hierarchical clustering tree of the matrix.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   C --  Coclassification matrix or other similarity matrix.
%
%
% Output Arguments
%__________________________________________________________________________
%
%   s -- Permutation of the nodes such that C(s,s) has large entries 
%        concentrated near the diagonal
%
% See also treeSort, drawHierarchy, consensusPlot

% Version: 1.1-alpha1
% Date: Tue Jan 16 18:15:01 EST 2018
% Author: Lucas Jeub
% Email: ljeub@iu.edu

% shortcut for single element input
if numel(C)==1
    s=1;
else
    
    % ensure symmetry
    C=(C+C')/2;
    
    % convert to distance
    C=full(max(C(:))-C);
    C=C-diag(diag(C));
    
    % hierarchical clustering
    Z=linkage(squareform(C),'average');
    
    % optimal order
    s=optimalleaforder(Z,C);
end
end
