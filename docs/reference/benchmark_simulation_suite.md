# Benchmark Simulation Methods Across Multiple Single-Cell Datasets

A unified, high-performance, and beginner-friendly benchmarking runner
designed for multi-dataset and multi-simulator studies (e.g., 16
datasets x 6 simulators). Automatically manages memory allocation with
sparse matrix coercions, parallel execution across datasets, progress
tracking, and produces both granular and cross-dataset aggregated
benchmark tables ready for visualization.

## Usage

``` r
benchmark_simulation_suite(
  datasets,
  simulators,
  compute_bivariate = FALSE,
  threads = 1,
  auto_sparse = TRUE,
  run_gc = TRUE,
  save_path = NULL,
  verbose = TRUE
)
```

## Arguments

- datasets:

  A named list of empirical reference datasets. Each item can be a
  single-cell count matrix (features x cells), a path to an RDS/CSV
  matrix file, or a list with `ref`, optional `cell_types`, and optional
  `batch_info`.

- simulators:

  A named list of simulation methods. Each element must be a named list
  of simulated count matrices (or matrix file paths) corresponding to
  the datasets in `datasets`. Example:
  `list(Splatter = list(D1 = sim_splat_d1, D2 = sim_splat_d2), scDesign3 = ...)`.

- compute_bivariate:

  Logical. Whether to compute 2D bivariate tests. Default `FALSE` for
  fast high-throughput execution.

- threads:

  Number of threads / parallel workers. Default is `1`.

- auto_sparse:

  Logical. Automatically coerce dense count matrices to `dgCMatrix` to
  minimize RAM usage. Default `TRUE`.

- run_gc:

  Logical. Run garbage collection
  ([`gc()`](https://rdrr.io/r/base/gc.html)) after each dataset
  evaluation to prevent memory accumulation on HPC nodes. Default
  `TRUE`.

- save_path:

  Optional file path (ending in `.rds`) to directly persist the
  consolidated benchmark bundle for loading into Shiny Studio.

- verbose:

  Logical. Print informative progress logging. Default `TRUE`.

## Value

An object of class `scSimEval_multidataset_benchmark` containing:

- `granular_benchmark_table`: Tidy data.frame containing all evaluated
  metrics across every dataset and simulator.

- `benchmark_summary_table`: Tidy data.frame containing cross-dataset
  averaged fidelity scores for each simulator (ready for
  [`plot_benchmark_bubble_matrix()`](https://kabilanbio.github.io/scSimEval/reference/plot_benchmark_bubble_matrix.md)).

- `dataset_summary_table`: Summary of reference dataset dimensions,
  sparsity, and library properties.

- `methods`: Vector of evaluated simulator names.

- `datasets`: Vector of evaluated dataset names.

- `cross_dataset_stats`: Summary statistics (mean, SD, SE) for each
  metric across datasets.

## Examples

``` r
# \donttest{
data(example_scrna, package = "scSimEval")
refs <- list("CohortA" = example_scrna$ref)
sims <- list("Splatter" = list("CohortA" = example_scrna$sim))
res <- benchmark_simulation_suite(refs, sims, verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
print(res)
#> =====================================================================
#>   scSimEval Multi-Dataset Benchmark Results
#> =====================================================================
#> Datasets Evaluated     : 1 (CohortA)
#> Simulators Evaluated   : 1 (Splatter)
#> Total Evaluated Scores : 166 across all datasets & metrics
#> Cross-Dataset Average  : 36 averaged metric entries in summary table
#> =====================================================================
# }
```
