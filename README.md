# HierarchicalConsensus


This code implements the hierarchical consensus clustering method introduced in

**Multiresolution Consensus Clustering in Networks**  
Lucas G. S. Jeub, Olaf Sporns, and Santo Fortunato  
[Scientific Reports 8, Article number: 3259 (2018)](https://doi.org/10.1038/s41598-018-21352-7)  
[arXiv:1710.02249](https://arxiv.org/abs/1710.02249)

Please cite this paper and include a link to https://github.com/LJeub/HierarchicalConsensus if you include results based on this code in an academic publication.


## Installation 

Download the latest version of the code [here](https://github.com/LJeub/HierarchicalConsensus/releases/latest) and add it to your MATLAB path. By default the code relies on [GenLouvain](http://netwiki.amath.unc.edu/GenLouvain/GenLouvain) to identify community structure. To use the code with default options make sure GenLouvain is installed on your MATLAB path (more information is available on the [GenLouvain GitHub page](https://github.com/GenLouvain/GenLouvain)). 


## Usage

This section gives a brief overview of the different functions included in this package. For more details use `help('function')` or `doc('function')`. The explanations below assume that we are given a network with adjacency matrix `A` as either a full or sparse matrix.

#### Generate initial ensemble

To use the hierarchical consensus clustering algorithm we first need to generate an ensemble of input partitions `S`, where `S(i,t)` gives the community membership of node `i` in partition `t`.  We use `np` for the number of partitions in the ensemble. The package includes two functions 

```Matlab
S = eventSamples(A, np)

S = exponentialSamples(A, np)
```

that compute ensembles based on multiresolution modularity. `eventSamples` implements the event sampling procedure and is usually recommended.
To generate an ensemble of input partitions fixing `gamma=1`, use

```Matlab
S = fixedResSamples(A, np)
```

For a different fixed resolution, use

```Matlab
S = fixedResSamples(A, np, 'Gamma', gamma)
```

Alternatively, one could use any other method to compute the initial ensemble. 

#### Hierarchical Consensus Community Structure

Given an initial ensemble `S` identify consensus community structure at significance level `alpha=0.05`:

```Matlab
[Sc, Tree] = hierarchicalConsensus(S);
```

For a different significance level `alpha` use

```Matlab
[Sc, Tree] = hierarchicalConsensus(S, alpha);
```

The `Tree` returned by `hierarchicalConsensus` is weighted, where the weight indicates the strength of the new cluster. This weight is used for visualizing the hierarchy and to extract partitions that are representative of different scales. By default the weight is the mean value of the coclassification matrix restricted to the nodes within the cluster. The way the weight is computed can be changed by using the options `'SimilarityFunction'` and `'SimilarityType'`. The weight has no influence on the way clusters are merged and it is possible that the weights are inconsistent with the hierarchical structure and the function issues a warning in that case. For more details and other optional parameters see `help('hierarchicalConsensus')`.

To change the weights for an existing hierarchy (e.g., because the existing weights are inconsistent), use
```Matlab
C = coclassificationMatrix(S);
[Tree, isConsistent] = dendrogramSimilarity(C, Sc, Tree, ...)
```

which accepts the same `'SimilarityFunction'` and `'SimilarityType'` options as `hierarchicalConsensus`. The second output argument `isConsistent` indicates if the weights are consistent with the hierarchy.

To compute all cuts of the consensus hierarchy that result from cutting the dendrogram use

```Matlab
[Sall, thresholds] = allPartitions(Sc, Tree)
```

The partition `Sall(:, i)` is the result of cutting the dendrogram at any threshold `t` such that `thresholds(i-1) < t <= thresholds(i)`. 


#### Visualize Results

Plot consensus hierarchy:

```Matlab
drawHierarchy(Sc, Tree)
```

Also show the coclassification matrix

```Matlab
C = coclassificationMatrix(S);
consensusPlot(C, Sc, Tree)
```

#### Generate Hierarchical Benchmark Networks

Choose the fraction of edges assigned to each level of the hierarchy by specifying `p`, where `p(i)` is the fraction of edges assigned to the `i-1`th level of the hierarchy (`p(1)` is the fraction of random edges not constrained by any community structure). Note that `sum(p)==1`.

Generate adjacency matrix `A` and ground truth partitions `Sgtruth`:

```Matlab
[A, Sgtruth] = hierarchicalBenchmark(n, p)
```

where `n` is the number of nodes in the network. For more customization options see `help('hierarchicalBenchmark')`.


## Example

This code produces a figure that is similar to Fig. 3a of the paper. 

First, generate the benchmark network:

```Matlab
[A, Sgtruth] = hierarchicalBenchmark(1000, [0.2,0.2,0.6]);
```

Next, compute the initial ensemble with 1000 partitions:

```Matlab
S = eventSamples(A, 1000);
```

Find hierarchical consensus community structure at significance level `alpha=0.05`:

```Matlab
[Sc, Tree] = hierarchicalConsensus(S);
```

Plot the results:

```Matlab
C = coclassificationMatrix(S);

consensusPlot(C, Sc, Tree, 'GroundTruth', Sgtruth)
```
