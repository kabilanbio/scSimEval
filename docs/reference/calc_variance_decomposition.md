# Variance Component Decomposition (BASiCS & muscat)

Decomposes gene expression variance into biological cell-type
heterogeneity, donor / biological replicate variation, and residual
technical noise.

## Usage

``` r
calc_variance_decomposition(
  counts,
  cell_metadata,
  donor_col = "donor",
  celltype_col = "cell_type",
  n_genes = 200
)
```

## Arguments

- counts:

  Matrix or data.frame (genes x cells) or SingleCellExperiment.

- cell_metadata:

  Data.frame containing donor and cell-type columns.

- donor_col:

  Character, column name for donor/replicate (default "donor").

- celltype_col:

  Character, column name for cell identity (default "cell_type").

- n_genes:

  Number of top variable genes to analyze (default 200).

## Value

A list of mean percentage of variance explained by cell-type, donor, and
residual noise.
