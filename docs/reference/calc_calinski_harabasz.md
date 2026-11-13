# Calculate Calinski-Harabasz Index (CH)

Calculate Calinski-Harabasz Index (CH)

## Usage

``` r
calc_calinski_harabasz(data, cluster_labels)
```

## Arguments

- data:

  Matrix with features in rows, cells in columns.

- cluster_labels:

  Cluster assignments.

## Value

Calinski-Harabasz index.

## Examples

``` r
data <- matrix(stats::rnorm(200), 20, 10)
cl <- factor(rep(c("A", "B"), each = 5))
dist_mat <- stats::dist(data)
calc_calinski_harabasz(data, cl)
#> [1] 0.8093449
```
