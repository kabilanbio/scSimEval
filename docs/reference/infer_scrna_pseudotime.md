# Automatically Infer Pseudotime Trajectory from scRNA-seq Counts

Estimates cell differentiation pseudotime directly from single-cell
expression counts using diffusion/principal curve projection along the
first principal component of top highly variable genes.

## Usage

``` r
infer_scrna_pseudotime(data, n_top = 500)
```

## Arguments

- data:

  Expression or count matrix (genes x cells).

- n_top:

  Number of top variable genes to use (default 500).

## Value

Numeric vector of pseudotime values in \[0, 1\] for each cell.
