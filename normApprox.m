function P=normApprox(alpha,mu)
% normApprox Compute normal approximate null model
%
% Syntax
%__________________________________________________________________________
%
%   P=normApprox(alpha,mu)
%
%
% Description
%__________________________________________________________________________
%
%   P=normApprox(alpha,mu) computes a normal approximation to the
%       alpha-confidence level of the Poisson-Binomial distribution with
%       parameters mu
%
%
% Input Arguments
%__________________________________________________________________________
%
%   alpha -- Confidence level to compute
%            0 <= alpha <= 1
%
%   mu -- Parameters for the Poisson-Binomial distribution. If
%       'size(mu,1)==1', the Poisson-Binomial distribution is assumed to be
%       the same for each node, otherwise each row of 'mu' specifies the
%       parameters for the Poisson-Binomial distribution for the
%       corresponding node.
%
%
% Output Arguments
%__________________________________________________________________________
%
%   P -- Null model at significance level 'alpha'. This is an 'NxN' matrix,
%        where 'size(mu,1)==N'. 
%
%
% Notes
%__________________________________________________________________________
%
% When the input 'size(mu,1)==1', P will be a scalar. This is the case for 
% the permutation model. If 'size(mu,1)>1', the function computes the
% confidence level for each row of mu. The output 'P(i,j)' is the minimum
% of the confidence level between rows 'i' and 'j'. The normal
% approximation takes into account that there is no randomness when
% 'mu==1' and the output is clipped to [0,1].
%
% See Also sampleApprox, localPermModel, permModel, hierarchicalConsensus

% Version: 1.0
% Date: Thu Oct  5 17:09:19 EDT 2017
% Author: Lucas Jeub
% Email: ljeub@iu.edu

[N,L]=size(mu);
p=zeros(N,1);
for i=1:N
    c=mu(i,:)==1;
    k_low=sum(c);
    
    sigma=sqrt((sum(mu(i,:)-mu(i,:).^2)));
    mui=sum(mu(i,~c));
    if sigma==0
        p(i)=(k_low+mui)/L;
    else
        p(i)=min(max(0,(k_low+norminv(alpha,mui,sigma))/L),1);
    end
    
end
P=zeros(N,N);
for i=1:N
    for j=1:N
        P(i,j)=min(p(i),p(j));
    end
end
end
