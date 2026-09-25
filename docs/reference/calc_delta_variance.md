# Calculate Delta Variance (Pseudoreplication Bias Metric)

Measures the difference in gene expression variance between true
biological replicates and cell-shuffled pseudoreplicates (Squair et al.,
Nature Communications 2021). Positive delta variance indicates that
biological replicates have higher inter-individual variability than
randomly pooled pseudoreplicates.

## Usage

``` r
calc_delta_variance(
  counts,
  replicates,
  conditions,
  cell_types = NULL,
  n_perm = 3
)
```

## Arguments

- counts:

  Count matrix (features x cells).

- replicates:

  Factor of biological replicate / batch assignments for each cell.

- conditions:

  Factor of condition or cell-type labels.

- cell_types:

  Optional cell type vector to evaluate per cell type or overall.

- n_perm:

  Number of permutation shuffles to compute mean pseudoreplicate
  variance. Default is 3.

## Value

A named list:

- delta_variance:

  Numeric vector of gene-level delta variance values

- mean_delta_variance:

  Mean delta variance across all genes

- median_delta_variance:

  Median delta variance across all genes
