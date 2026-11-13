# Full Clustering Performance Evaluation

Evaluates both unsupervised cluster separation (ASW, DB, Neighborhood
Purity, CDI) and supervised concordance against ground truth labels
(ARI, NMI).

## Usage

``` r
evaluate_clustering_metrics(
  data,
  cluster_info = NULL,
  pred_clusters = NULL,
  dist_mat = NULL,
  cell_types = NULL,
  ref_data = NULL,
  ref_cell_types = NULL,
  batch_info = NULL,
  threads = 1L
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

- ref_cell_types:

  Optional reference cell type labels (if ref_data has different cells
  than data).

- batch_info:

  Optional batch labels per cell passed to
  [`calc_cdi`](https://kabilanbio.github.io/scSimEval/reference/calc_cdi.md).

- threads:

  Integer number of parallel CPU threads passed to
  [`calc_cdi`](https://kabilanbio.github.io/scSimEval/reference/calc_cdi.md)
  (default `1L`).

## Value

A named list of all clustering metrics.

## Examples

``` r
data <- matrix(stats::rpois(200, lambda = 3), 20, 10)
cl <- factor(rep(c("A", "B"), each = 5))
dist_mat <- stats::dist(t(data))
evaluate_clustering_metrics(data, cl)
#> $silhouette_sim
#> [1] 0.07265493
#> 
#> $silhouette
#> [1] 0.07265493
#> 
#> $davies_bouldin
#> [1] 2.682595
#> 
#> $davies_bouldin_sim
#> [1] 2.682595
#> 
#> $neighborhood_purity
#> [1] 0.46
#> 
#> $cdi
#> [1] 4.292815
#> 
#> $CDI_AIC
#> [1] 4.292815
#> 
#> $CDI_BIC
#> [1] 4.413849
#> 
#> $ari
#> [1] 0
#> 
#> $ARI
#> [1] 0
#> 
#> $nmi
#>         1 
#> 0.1470822 
#> 
#> $NMI
#>         1 
#> 0.1470822 
#> 
```
