# Full Clustering Performance Evaluation

Evaluates both unsupervised cluster separation (ASW, Dunn, Connectivity,
DB, CH, Neighborhood Purity, CDI) and supervised concordance against
ground truth labels (Clustering Accuracy ACC, Hungarian
F1/Precision/Recall, ARI, NMI, AMI, FMI, Homogeneity, Completeness,
V-measure).

## Usage

``` r
evaluate_clustering_metrics(
  data,
  cluster_info = NULL,
  pred_clusters = NULL,
  dist_mat = NULL,
  cell_types = NULL,
  ref_data = NULL
)
```

## Arguments

- data:

  Count or normalized matrix (features x cells).

- cluster_info:

  Cluster or cell type labels.

- pred_clusters:

  Optional predicted cluster labels (if different from ground truth). If
  NULL, k-means is automatically run.

- dist_mat:

  Optional precomputed distance matrix.

- cell_types:

  Optional alias for cluster_info.

- ref_data:

  Optional reference data to compute reference clustering quality.

## Value

A named list of all clustering metrics.
