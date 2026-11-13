# Construct Pseudobulk Count Matrices from Single-Cell Data (Libra Engine)

Converts a single-cell count matrix (features x cells) for scRNA-seq or
scATAC-seq to sample-level pseudobulk matrices by aggregating (summing)
single-cell counts within biological replicates and experimental
conditions per cell type (Squair et al., Nature Communications 2021 /
Libra).

## Usage

``` r
to_pseudobulk(
  input,
  cell_types,
  replicates,
  conditions = NULL,
  binarize = FALSE,
  min_cells = 3,
  min_reps = 2,
  min_features = 0
)
```

## Arguments

- input:

  Single-cell count matrix (features x cells), dgCMatrix,
  SingleCellExperiment, or Seurat object.

- cell_types:

  Factor or vector specifying cell type annotations for each cell.

- replicates:

  Factor or vector specifying biological replicate / batch annotations
  for each cell. Required.

- conditions:

  Optional factor or vector specifying experimental conditions or
  comparison labels.

- binarize:

  Logical. If TRUE, binarizes single-cell counts (cells \> 0 set to 1)
  prior to aggregation, producing accessibility-count pseudobulks
  (standard for scATAC-seq in Libra). Default is FALSE.

- min_cells:

  Minimum number of cells required per cell type to retain. Default is
  3.

- min_reps:

  Minimum number of replicates required to retain. Default is 2.

- min_features:

  Minimum non-zero count threshold to retain a feature. Default is 0.

## Value

A named list of pseudobulk count matrices (features x
sample_replicates), one for each retained cell type.

## Examples

``` r
data(example_scatac, package = "scSimEval")
rep_id <- factor(rep(c("Rep1", "Rep2"), length.out = ncol(example_scatac$ref)))
ct_id <- factor(rep(c("TypeA", "TypeB"), each = ncol(example_scatac$ref) / 2))
pb <- to_pseudobulk(example_scatac$ref, cell_types = ct_id, replicates = rep_id)
length(pb)
#> [1] 2
```
