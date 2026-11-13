# Calculate Davies-Bouldin Index (DB)

Calculate Davies-Bouldin Index (DB)

## Usage

``` r
calc_davies_bouldin(data, cluster_labels)
```

## Arguments

- data:

  Matrix with features in rows, cells in columns.

- cluster_labels:

  Cluster assignments.

## Value

Davies-Bouldin index (lower is better).

## Examples

``` r
data <- matrix(stats::rnorm(200), 20, 10)
cl <- factor(rep(c("A", "B"), each = 5))
dist_mat <- stats::dist(data)
calc_davies_bouldin(data, cl)
#> [1] 3.091148
```
