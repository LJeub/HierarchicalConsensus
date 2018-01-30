function mu=localPermModel(S)
% localPermModel Compute co-classification probabilities using local Permutation Model
%
% Syntax
%__________________________________________________________________________
%
%   mu=localPermModel(S)
%
%
% Description
%__________________________________________________________________________
%
%   mu=localPermModel(S) computes co-classification probabilities for each
%       node and partition based on a local permutation model.
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
%   mu -- Matrix of probabilities, where 'mu(i,j)' is the probability for
%         node 'i' to be co-assigned with any other node in partition 'j'
%         given that the community assignment of 'i' remains fixed.
%
% See Also hierarchicalConsensus, permModel, normApprox, sampleApprox

% Version: 1.1
% Date: Tue 30 Jan 2018 18:22:38 EST
% Author: Lucas Jeub
% Email: ljeub@iu.edu


[N,L]=size(S);
C=max(S(:));
mu=zeros(N,L);
if N>1
    sizes=zeros(C,L);
    for i=1:L
        G=sparse(1:N,S(:,i),1,N,C);
        sizes(:,i)=(sum(G,1));
    end
    
    for i=1:L
        mu(:,i)=(sizes(S(:,i),i)-1)./(N-1);
    end
end
end
