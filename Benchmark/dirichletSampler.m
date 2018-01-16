function weights=dirichletSampler(theta,nc)
% dirichletSampler Sample weights from a symmetric Dirichlet distribution
%
% Syntax
%__________________________________________________________________________
%
%   weights=dirichletSampler(theta,nc)
%
%
% Description
%__________________________________________________________________________
%
%   weights=dirichletSampler(theta,nc) samples a discrete probability
%       distribution on 'nc' elements from a symmetric Dirichlet 
%       distribution with shape parameter 'theta'.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   theta -- shape parameter for Dirichlet distribution
%
%   nc -- number of weights to sample
%
%
% Output Arguments
%__________________________________________________________________________
%
%   weights -- sampled weigths (i.e. probabilities of a categorical
%       distribution)
%
%
% See also hierarchicalBenchmark

% Version: 1.1-alpha1
% Date: Tue Jan 16 18:15:01 EST 2018
% Author: Lucas Jeub
% Email: ljeub@iu.edu

weights=gamrnd(theta,1,nc,1);
weights=weights./sum(weights);

end

