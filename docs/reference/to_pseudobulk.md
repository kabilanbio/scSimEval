# Convert Single-Cell Counts to Replicate Pseudobulks

Converts a single-cell expression matrix into pseudobulk matrices per
cell type by summarizing counts across biological replicates and
experimental conditions (adapted from Squair et al. Libra).

## Usage

``` r
to_pseudobulk(counts, cell_types, replicates, conditions, min_cells = 3)
```

## Arguments

- counts:

  Single-cell count matrix (features x cells).

- cell_types:

  Factor or character vector of cell types.

- replicates:

  Factor or character vector of biological replicate IDs.

- conditions:

  Factor or character vector of experimental condition / treatment
  labels.

- min_cells:

  Minimum number of cells required per cell type to retain. Default is
  3.

## Value

Named list of pseudobulk count matrices (features x sample_replicates),
one per cell type.
