# HierarchicalConsensus


This code implements the hierarchical consensus clustering method introduced in

**Multiresolution Consensus Clustering in Networks**  
Lucas G. S. Jeub, Olaf Sporns, and Santo Fortunato  
[arXiv:xxxx.xxxx](http://arxiv.com/)

Please cite this paper and include a link to https://github.com/LJeub/HierarchicalConsensus if you include results based on this code in an academic publication.


## Installation 

Download the latest version of the code [here](https://github.com/LJeub/HierarchicalConsensus/releases/latest) and add it to your MATLAB path. By default the code relies on [GenLouvain](http://netwiki.amath.unc.edu/GenLouvain/GenLouvain) to identify community structure. To use the code with default options make sure GenLouvain is installed on your MATLAB path (more information is available on the [GenLouvain GitHub page](https://github.com/GenLouvain/GenLouvain)). 


## Usage

This section gives a brief overview of the different functions included in this package. For more details use `help('function')` or `doc('function')`. The explanations below assume that we are given a network with adjacency matrix `A` as either a full or sparse matrix.

#### Generate initial ensemble

To use the hierarchical consensus clustering algorithm we first need to generate an ensemble of input partitions`S`, were `S(i,t)` gives the community membership of node `i` in partition `t`.  The package includes two functions 

	S = eventSamples(A, np)

	S = exponentialSamples(A, np)
	
that compute ensembles based on multiresolution modularity. `eventSamples` implements the event sampling procedure and is usually recommended. Alternatively, one could use any other method to compute the initial ensemble. Here `np` is the number of partitions in the ensemble. 

#### Hierarchical Consensus Community Structure

Given an initial ensemble `S` identify consensus community structure at significance level `alpha=0.05`:

	[Sc, Tree] = hierarchicalConsensus(S);

For a different significance level `alpha` use

	[Sc, Tree] = hierarchicalConsensus(S, alpha);
	
For more details other optional parameters see `help('hierarchicalConsensus')`.

Compute all cuts of the consensus hierarchy for further analysis:
	
	Sall = allPartitions(Sc, Tree)

#### Visualize Results

Plot consensus hierarchy:
	
	drawHierarchy(Sc, Tree)

Also show the coclassification matrix

	C = coclassificationMatrix(S);
	consensusPlot(C, Sc, Tree)
	
#### Generate Hierarchical Benchmark Networks

Choose the fraction of edges assigned to each level of the hierarchy by specifying `p`, where `p(i)` is the fraction of edges assigned to the `i+1`th level of the hierarchy (`p(1)` is the fraction of random edges not constraint by any community structure). Note that `sum(p)==1`.

Generate adjacency matrix `A` and ground truth partitions `Sgtruth`:

	[A, Sgtruth] = hierarchicalBenchmark(n, p)
	
where `n` is the number of nodes in the network. For more customization options see `help('hierarchicalBenchmark')`.


## Example

This code produces a figure that is similar to Fig. 3a of the paper. 

First, generate the benchmark network:

	[A, Sgtruth] = hierarchicalBenchmark(1000, [0.2,0.2,0.6]);
	

Next, compute the initial ensemble with 150 partitions:

	S = eventSamples(A, 150);
	
Find hierarchical consensus community structure at significance level `alpha=0.05`:

	[Sc, Tree] = hierarchicalConsensus(S);
	
Plot the results:
	
	C = coclassificationMatrix(S);
	
	consensusPlot(C, Sc, Tree, 'GroundTruth', Sgtruth)
