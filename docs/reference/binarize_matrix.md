# Binarize Single-Cell Count Matrix (scATAC-seq Feature Binarization)

Sets all non-zero counts in a single-cell feature matrix to 1 (features
x cells), converting quantitative read cut-counts into binary detection
status. Frequently applied to scATAC-seq peak accessibility data prior
to pseudobulking or single-cell modeling (Squair et al., 2021; Teo et
al., 2024 / Libra).

## Usage

``` r
binarize_matrix(mat)
```

## Arguments

- mat:

  Single-cell matrix, sparse dgCMatrix, SingleCellExperiment, or Seurat
  object.

## Value

A matrix of the same class with all positive entries set to 1.

## Examples

``` r
m <- matrix(c(0, 2, 1, 0, 5, 0), nrow = 2)
binarize_matrix(m)
#>      [,1] [,2] [,3]
#> [1,]    0    1    1
#> [2,]    1    0    0
```
