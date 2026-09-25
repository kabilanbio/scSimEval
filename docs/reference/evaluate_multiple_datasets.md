# Evaluate Multiple Single-Cell Datasets Simultaneously

Orchestrates simultaneous benchmarking across a collection of
single-cell datasets, including paired multiomics (e.g., scRNA-seq and
scATAC-seq from the same biological experiment) and unimodal datasets
(e.g., standalone scRNA-seq or scATAC-seq). Produces a unified,
consolidated summary table, comparative ranking matrix, and dataset
overview.

## Usage

``` r
evaluate_multiple_datasets(
  datasets,
  pair_by_prefix = TRUE,
  compute_bivariate = FALSE,
  threads = 1,
  verbose = TRUE
)
```

## Arguments

- datasets:

  A named list of datasets to evaluate. Can be provided in two formats:

  - **Format A (Grouped Multiomics & Unimodal):** A list of dataset
    specifications, e.g.:
    `list("Method 1" = list(ref = list(rna = ..., atac = ...), sim = list(rna = ..., atac = ...)), "Method 3" = list(ref = ..., sim = ..., modality = "scATAC-seq"))`

  - **Format B (Flat Named Datasets):** A list with named modality
    datasets, e.g.:
    `list("Method 1-scRNA-seq" = list(ref = ..., sim = ...), "Method 1-scATAC-seq" = list(ref = ..., sim = ...), "Method 3-scATAC-seq" = list(ref = ..., sim = ...), "Method 4-scRNA-seq" = list(ref = ..., sim = ...))`

- pair_by_prefix:

  Logical; if `TRUE` (default), datasets sharing a common root prefix
  (such as "Method 1-scRNA-seq" and "Method 1-scATAC-seq") are
  automatically paired as multiomics for joint cross-modal evaluation.
  If `FALSE`, each dataset entry is evaluated independently.

- compute_bivariate:

  Logical; whether to compute 2D bivariate tests (default FALSE for
  speed).

- threads:

  Number of CPU threads for parallel metric computation. Default is 1.

- verbose:

  Logical; whether to print progress messages. Default is TRUE.

## Value

An S3 object of class `scSimEval_consolidated` containing:

- `consolidated_summary_table`: Master tidy data.frame with columns
  `Dataset`, `Data_Name`, `Modality`, `Category`, `Property`, `Metric`,
  and `Value`.

- `consolidated_score_matrix`: Wide comparison matrix comparing all
  datasets across metrics.

- `dataset_overview`: Aggregate summary table with total metrics and
  mean scores per dataset.

- `dataset_results`: Named list containing detailed evaluation results
  for each individual dataset.
