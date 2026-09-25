# Extract Comprehensive Feature-Level (Gene / Peak) Properties

Computes mean abundance, standard deviation, variance, coefficient of
variation (CV), feature dropout (zero fraction), dispersion, Biological
Coefficient of Variation (BCV), gene-gene correlations, and outlier
proportions.

## Usage

``` r
extract_feature_properties(data, n_top_cor = 400, verbose = FALSE)
```

## Arguments

- data:

  Count matrix (features x cells).

- n_top_cor:

  Number of top HVGs to use for gene-gene correlation matrix (default
  400).

- verbose:

  Logical, whether to print progress messages.

## Value

A named list of feature-level summary vectors.
