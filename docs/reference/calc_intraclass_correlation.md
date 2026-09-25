# Calculate Intraclass Correlation Coefficient (ICC) (hierarchicell & rescueSim)

Computes the one-way random effects Intraclass Correlation Coefficient
(ICC) per gene across donors or biological replicates to quantify
hierarchical clustering.

## Usage

``` r
calc_intraclass_correlation(counts, donor_labels, n_genes = 200)
```

## Arguments

- counts:

  Matrix (genes x cells) or SingleCellExperiment.

- donor_labels:

  Vector or factor of donor / subject IDs.

- n_genes:

  Number of top variable genes to evaluate (default 200).

## Value

A list with median ICC, mean ICC, and vector of per-gene ICC values.
