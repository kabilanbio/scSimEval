# Calculate Clustering Accuracy (ACC)

Integrated from scCluBench (Xu et al., AAAI 2026). Aligns predicted
clusters to ground truth classes using the Hungarian bipartite matching
algorithm, then computes the fraction of correctly mapped cells.

## Usage

``` r
calc_clustering_accuracy(pred, truth)
```

## Arguments

- pred:

  Vector of predicted cluster labels.

- truth:

  Vector of ground truth cell type labels.

## Value

Clustering accuracy between 0 and 1.
