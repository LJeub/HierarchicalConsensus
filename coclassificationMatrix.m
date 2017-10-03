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

% Version:
% Date:
% Author:
% Email:

% identify unique rows
[N,L]=size(S);
[ur,~,ri]=unique(S,'rows');
Cr=zeros(size(ur,1));

% compute scores
for i=1:size(ur,1)
    Cr(i,:)=sum(bsxfun(@eq,ur(i,:),ur),2)/L;
end

% expand
G=sparse(1:N,ri,1);
C=G*Cr*G';

end