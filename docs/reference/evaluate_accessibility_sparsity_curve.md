# Evaluate Accessibility-Sparsity Curve Concordance (simATAC)

Compares the non-linear relationship between peak mean accessibility and
non-zero proportion (NZP) between reference and simulated scATAC-seq
datasets.

## Usage

``` r
evaluate_accessibility_sparsity_curve(ref_data, sim_data, poly_degree = 2)
```

## Arguments

- ref_data:

  Reference count matrix (peaks x cells).

- sim_data:

  Simulated count matrix (peaks x cells).

- poly_degree:

  Degree of polynomial (default: 2).

## Value

A list of reference and simulated curve parameters, absolute
discrepancies, and curve prediction RMSE.
