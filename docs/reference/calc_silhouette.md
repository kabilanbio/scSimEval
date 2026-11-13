# Calculate Average Silhouette Width (ASW)

Calculate Average Silhouette Width (ASW)

## Usage

``` r
calc_silhouette(dist_mat, cluster_labels)
```

## Arguments

- dist_mat:

  Distance matrix among cells or numeric matrix (features x cells).

- cluster_labels:

  Vector of cluster or cell-type labels.

## Value

Mean silhouette width.

## Examples

``` r
data <- matrix(stats::rnorm(200), 20, 10)
cl <- factor(rep(c("A", "B"), each = 10))
dist_mat <- stats::dist(data)
calc_silhouette(dist_mat, cl)
#> [1] -0.0127335
```
