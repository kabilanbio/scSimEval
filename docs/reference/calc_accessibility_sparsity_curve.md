# Fit Chromatin Accessibility-Sparsity Polynomial Curve (simATAC)

Fits a polynomial curve relating peak/bin mean accessibility to non-zero
cell proportion (NZP / detection frequency) as modeled in the simATAC
framework (Navidi et al., Genome Biology 2021): NZP = c0 + c1 \* mean +
c2 \* mean^2.

## Usage

``` r
calc_accessibility_sparsity_curve(data, poly_degree = 2)
```

## Arguments

- data:

  Count matrix (features/peaks x cells) or data.frame with 'mean' and
  'nzp'.

- poly_degree:

  Degree of polynomial (default: 2 for quadratic curve).

## Value

A list with estimated coefficients (c0, c1, c2), R-squared,
Spearman/Pearson correlations, and model fit summary.
