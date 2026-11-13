# Validate and Standardize Single-Cell Input Objects

Strictly enforces standard input object compliance in `scSimEval`
(Option B). Validates that input datasets are either
`SingleCellExperiment` or `Seurat` objects, extracting underlying count
matrices, cell annotations, batch labels, and dimensional reductions.

## Usage

``` r
validate_sc_object(obj, arg_name = "data", reduction_name = NULL)
```

## Arguments

- obj:

  An object to validate, strictly required to be a
  `SingleCellExperiment` or `Seurat` object.

- arg_name:

  Name of the argument being validated (for clear error messages).

- reduction_name:

  Optional name of precomputed dimensional reduction (e.g. `"PCA"`,
  `"integrated"`).

## Value

A list with standardized components:

- object:

  The original `SingleCellExperiment` or `Seurat` object.

- counts:

  Extracted raw count matrix (features x cells).

- cell_types:

  Extracted or parsed cell-type factor/vector, or NULL if absent.

- batch_info:

  Extracted or parsed batch factor/vector, or NULL if absent.

- embedding:

  Extracted dimensional reduction matrix (cells x dims), or NULL if
  absent.
