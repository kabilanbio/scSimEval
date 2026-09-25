# Evaluate Comprehensive Simulation Accuracy (Unimodal Omics Layer)

Computes distribution distances, cellular and feature-level properties,
zero-inflation curves, and manifold distances between reference and
simulated datasets.

## Usage

``` r
evaluate_simulation_accuracy(
  ref_data,
  sim_data,
  compute_bivariate = TRUE,
  threads = 1,
  cpu_time = NULL,
  memory_mb = NULL,
  system_time = NULL,
  elapsed_time = NULL,
  verbose = TRUE
)
```

## Arguments

- ref_data:

  Reference count matrix (features x cells).

- sim_data:

  Simulated count matrix (features x cells).

- compute_bivariate:

  Logical, whether to compute 2D bivariate tests. Default is TRUE.

- threads:

  Number of threads for parallel computation. Default is 1.

- cpu_time:

  Optional numeric value of CPU execution time in seconds.

- memory_mb:

  Optional numeric value of peak memory allocation in MB.

- system_time:

  Optional numeric value of system CPU time in seconds.

- elapsed_time:

  Optional numeric value of wall-clock elapsed time in seconds.

- verbose:

  Logical, whether to print execution messages. Default is TRUE.

## Value

A list containing tidy summary tables and detailed metric lists.
