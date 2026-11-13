# Validate and Standardize Multiomics Input Objects

Strictly enforces standard multiomics input compliance in `scSimEval`.
Validates that multiomics inputs conform to `Seurat` or
`SingleCellExperiment` structures. Accepts either:

- A named list of 2 `Seurat` or `SingleCellExperiment` objects (e.g.
  `list(rna = ..., atac = ...)`).

- A multi-modal `Seurat` object containing both RNA and ATAC assays.

- A multi-modal `SingleCellExperiment` object with an `altExp` (e.g.
  ATAC).

## Usage

``` r
validate_multiomics_input(obj, arg_name = "multi_data", reduction_name = NULL)
```

## Arguments

- obj:

  Multiomics dataset object to validate.

- arg_name:

  Name of argument.

- reduction_name:

  Optional reduction name to extract.

## Value

A list with standardized components: `m1`, `m2`, `cell_types`,
`batch_info`, `embedding`.
