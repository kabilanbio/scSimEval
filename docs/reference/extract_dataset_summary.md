# Extract Dataset Dimensional Properties and Summary Statistics

Computes key dimensional and biological properties from single-cell or
multiomics count matrices, including cell counts, feature counts,
sparsity percentage, biological group/cell type counts, batch counts,
median library size, median detected features per cell, and mean
expression level.

## Usage

``` r
extract_dataset_summary(
  mat,
  role = "Reference",
  method_name = "Dataset",
  modality = "scRNA-seq",
  cell_types = NULL,
  batch_info = NULL
)
```

## Arguments

- mat:

  A count matrix (genes/features x cells), dgCMatrix,
  SingleCellExperiment, or Seurat object.

- role:

  Character string describing the role of the dataset (e.g.
  `"Biological Reference"` or `"Simulated"`). Default is `"Reference"`.

- method_name:

  Character string indicating the method or dataset name (e.g.
  `"Empirical Reference"`, `"Splatter"`). Default is `"Dataset"`.

- modality:

  Character string describing the molecular modality (e.g.
  `"scRNA-seq"`, `"scATAC-seq"`). Default is `"scRNA-seq"`.

- cell_types:

  Optional factor or character vector of cell type annotations.

- batch_info:

  Optional factor or character vector of batch annotations.

## Value

A data.frame with 1 row summarizing the dataset properties.

## Examples

``` r
data(example_scrna, package = "scSimEval")
extract_dataset_summary(example_scrna$ref, role = "Reference", method_name = "Empirical")
#>   Dataset / Simulator      Role  Modality Cells (N) Features (P) Sparsity
#> 1           Empirical Reference scRNA-seq        80           60   17.92%
#>   Cell Types (Groups) Batches Median Lib Size Median Detected Features
#> 1                   2       2           257.0                       49
#>   Mean Expression
#> 1           4.272
```
