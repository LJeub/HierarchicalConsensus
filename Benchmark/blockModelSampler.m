function A=blockModelSampler(W,S,k)
% blockModelSampler Sample adjacency matrix for degree-corrected block-model
%
% Syntax
%__________________________________________________________________________
%
%   A=blockModelSampler(W,S,k)
%
%
% Description
%__________________________________________________________________________
%
%   A=blockModelSampler(W,S,k) samples an adjacency matrix based on a
%       degree-corrected stochastic block-model using binomial sampling for
%       dense blocks and rejection sampling to avoid multi-edges for sparse
%       blocks.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   W -- block-matrix giving expected number of edges between blocks
%   S -- block assignment vector
%   k -- node weights (will be normalized by total community weights)
%
%
% Output Arguments
%__________________________________________________________________________
%       
%   A -- sampled adjacency matrix
%
%
% See also hierarchicalBenchmark

% Version: 1.1-alpha1
% Date: Tue Jan 16 18:15:01 EST 2018
% Author: Lucas Jeub
% Email: ljeub@iu.edu

max_reject=100;
n_nodes=numel(S);
n_groups=max(S);
if n_nodes~=numel(k)
    error('dimension mismatch')
end
if n_groups~=size(W,1)||n_groups~=size(W,2)
    error('block-number mismatch')
end

G=sparse(1:n_nodes,S,true,n_nodes,n_groups);
group_sizes=sum(G,1);
k_g=arrayfun(@(i) sum(k(G(:,i))), 1:n_groups);
mm=sum(k);
pos=0;
neighbors=cell(n_nodes,1);
for i=1:n_nodes
    neighbors{i}=i;
end
sigma=cell(n_groups,1);
for i=1:n_groups
    sigma{i}=cumsum(k(G(:,i)))./k_g(i);
end

for g1=1:n_groups
    pos_blocks=g1:n_groups;
    for g2=pos_blocks(W(g1:n_groups,g1)>0)
        m=poissrnd(full(W(g1,g2)));
        nodes_g1=find(G(:,g1))';
        nodes_g2=find(G(:,g2))';
        if g1==g2
            dense=2*m>group_sizes(g1)*(group_sizes(g2)-1);
            m=m/2;
        else
            dense=2*m>group_sizes(g1)*group_sizes(g2);
        end
        pos=pos+m;
        
        if dense
            sig1=k(G(:,g1))./k_g(g1);
            sig2=k(G(:,g2))./k_g(g2);
            ng1=group_sizes(g1);
            ng2=group_sizes(g2);
            P=sig1*W(g1,g2)*sig2';
            if g1==g2
                for it=1:ng1
                    nid=find(P(it,it+1:end)>rand(1,ng1-it))+it;
                    neighbors{nodes_g1(it)}=...
                        [neighbors{nodes_g1(it)},nodes_g2(nid)];
                    for it2=1:length(nid)
                        neighbors{nodes_g2(nid(it2))}=...
                            [neighbors{nodes_g2(nid(it2))},nodes_g1(it)];
                    end
                end
            else
                for it=1:ng1
                    nid=find(P(it,:)>rand(1,ng2));
                    neighbors{nodes_g1(it)}=...
                        [neighbors{nodes_g1(it)},nodes_g2(nid)];
                    for it2=1:length(nid)
                        neighbors{nodes_g2(nid(it2))}=...
                            [neighbors{nodes_g2(nid(it2))},nodes_g1(it)];
                    end
                end
            end
        else
            for e=1:m
                isneighbor=true;
                reject_count=0;
                while isneighbor&&reject_count<=max_reject
                    decide=rand(2,1);
                    n1=nodes_g1(find(sigma{g1}>decide(1),1));
                    n2=nodes_g2(find(sigma{g2}>decide(2),1));
                    isneighbor=any(n1==neighbors{n2}(:));
                    reject_count=reject_count+1;
                end
                if reject_count>max_reject
                    warning('stopping sampling of edges for current block')
                    break;
                end
                neighbors{n1}=[neighbors{n1},n2];
                neighbors{n2}=[neighbors{n2},n1];
            end
        end
    end
end

indrow=zeros(2*pos,1);
indcol=zeros(2*pos,1);
pos=0;
for i=1:n_nodes
    indrow(pos+1:pos+length(neighbors{i})-1)=i;
    indcol(pos+1:pos+length(neighbors{i})-1)=neighbors{i}(2:end);
    pos=pos+length(neighbors{i})-1;
end
indrow=indrow(1:pos);
indcol=indcol(1:pos);
A=sparse(indrow,indcol,1,n_nodes,n_nodes);
end

            

