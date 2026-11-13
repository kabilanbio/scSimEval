# Calculate Dunn Index

Calculate Dunn Index

## Usage

``` r
calc_dunn(dist_mat, cluster_labels)
```

## Arguments

- dist_mat:

  Distance matrix or numeric matrix.

- cluster_labels:

  Cluster assignments.

## Value

Dunn index (ratio of smallest inter-cluster to largest intra-cluster
distance).

## Examples

``` r
data <- matrix(stats::rnorm(200), 20, 10)
cl <- factor(rep(c("A", "B"), each = 10))
dist_mat <- stats::dist(data)
calc_dunn(dist_mat, cl)
#> [1] 0.2504987
```
