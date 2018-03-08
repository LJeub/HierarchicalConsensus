function C=coclassificationMatrix(S)
% coclassificationMatrix Compute coclassification frequency between nodes
%
% Syntax
%__________________________________________________________________________
%
%   C=coclassificationMatrix(S)
%
%
% Description
%__________________________________________________________________________
%
%   C=coclassificationMatrix(S) computes the frequency each pair of nodes
%       appeared in the same cluster in an ensemble of partitions 'S'
%
%
% Input Arguments:
%__________________________________________________________________________
%
%   S -- Ensemble of partitions given as a matrix where each column
%        corresponds to a partition.
%
%
% Output Arguments:
%
%   C -- Coclassification matrix. C(i,j) encodes the frequency nodes 'i'
%        and 'j' appeared together in a cluster in S.
%
% See also hierarchicalConsensus

% Version: 1.1.1
% Date: Thu  8 Mar 2018 15:34:46 CET
% Author: Lucas Jeub
% Email: ljeub@iu.edu

[N,L]=size(S);
% reduce to unique rows
[S,~,ri]=unique(S,'rows');
E=sparse(1:N,ri,1);
% construct node-to-cluster adjacency matrix
row=repmat(1:size(S,1),1,L);
for i=2:L
    S(:,i)=S(:,i)+max(S(:,i-1));
end
G=sparse(row,S(:),1);
% compute coclassification matrix
C=E*(G*G')*E';
C=full(C)./L;
end
