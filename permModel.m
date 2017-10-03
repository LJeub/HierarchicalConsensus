function mu=permModel(S)
% permModel Compute co-classification probabilities using Permutation Model
%
% Syntax
%__________________________________________________________________________
%
%   mu=permModel(S)
%
%
% Description
%__________________________________________________________________________
%
%   mu=permModel(S) computes co-classification probabilities for each
%       node and partition based on the permutation model.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   S -- Ensemble of input partitions given as a matrix where each column
%   encodes a partition
%
%
% Output Arguments
%__________________________________________________________________________
%
%   mu -- Vector of probabilities, where 'mu(j)' is the probability for
%         a pair of nodes to be co-assigned in partition 'j' under the
%         permutation model.
%
% See Also hierarchicalConsensus, localPermModel, normApprox, sampleApprox

% Version:
% Date:
% Author:
% Email:

[N,L]=size(S);
C=max(S(:));
if N==1
    mu=zeros(N,L);
else
    sizes=zeros(C,L);
    for i=1:L
        G=sparse(1:N,S(:,i),1,N,C);
        sizes(:,i)=(sum(G,1))';
    end
    mu=sum((sizes.*(sizes-1))./(N*(N-1)));
end

end