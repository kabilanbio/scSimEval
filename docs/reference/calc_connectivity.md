# Calculate Cluster Connectivity

Calculate Cluster Connectivity

## Usage

``` r
calc_connectivity(dist_mat, cluster_labels)
```

## Arguments

- dist_mat:

  Distance matrix or numeric matrix.

- cluster_labels:

  Cluster assignments.

## Value

Connectivity metric.

## Examples

``` r
data <- matrix(stats::rnorm(200), 20, 10)
cl <- factor(rep(c("A", "B"), each = 10))
dist_mat <- stats::dist(data)
calc_connectivity(dist_mat, cl)
#> [1] 33.95556
```
