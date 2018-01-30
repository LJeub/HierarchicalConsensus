function P=sampleApprox(alpha,mu,n_samples)
% sampleApprox Compute sample-approximated null model
%
% Syntax
%__________________________________________________________________________
%
%   P=sampleApprox(alpha,mu)
%
%   P=sampleApprox(alpha,mu,n_samples)
%
%
% Description
%__________________________________________________________________________
%
%   P=sampleApprox(alpha,mu) computes a sampled approximation to the
%       alpha-confidence level of the Poisson-Binomial distribution with
%       parameters mu using 10000 samples for each partition
%
%   P=sampleApprox(alpha,mu,n_samples) uses 'n_samples' number of samples
%       for each partition.
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
%   n_samples -- Number of samples to use to approximate the distribution
%       (default: 10000)
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
% of the confidence level between rows 'i' and 'j'.
%
% See Also normApprox, localPermModel, permModel, hierarchicalConsensus

% Version: 1.1
% Date: Tue 30 Jan 2018 18:22:38 EST
% Author: Lucas Jeub
% Email: ljeub@iu.edu

if nargin<3
    n_samples=10000;
end
[N,L]=size(mu);

p=zeros(N,1);
r=rand(n_samples,L);
for i=1:N
    if sum(mu(i,:)-mu(i,:).^2)==0 % variance is 0
        p(i)=sum(mu(i,:))/L;
    else
        p(i)=quantile(sum(bsxfun(@lt,r,mu(i,:)),2),alpha)/L;
    end
end
P=zeros(N,N);
for i=1:N
    for j=1:N
        P(i,j)=min(p(i),p(j));
    end
end

end
