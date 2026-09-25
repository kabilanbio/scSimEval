# Calculate Correlation of Geodesic Pseudotime Distances

Evaluates whether the relative cell-to-cell ordering along a
differentiation trajectory is preserved between reference and simulated
data. Accepts either precomputed pseudotime vectors OR scRNA-seq
expression matrices (which will be automatically inferred via PCA
trajectory).

## Usage

``` r
calc_pseudotime_correlation(
  ref_pseudotime,
  sim_pseudotime,
  method = c("spearman", "pearson")
)
```

## Arguments

- ref_pseudotime:

  Numeric vector of reference pseudotimes OR reference count matrix
  (genes x cells).

- sim_pseudotime:

  Numeric vector of simulated pseudotimes OR simulated count matrix
  (genes x cells).

- method:

  Correlation method: "spearman" (default) or "pearson".

## Value

Correlation coefficient between reference and simulated pseudotime
trajectories.
