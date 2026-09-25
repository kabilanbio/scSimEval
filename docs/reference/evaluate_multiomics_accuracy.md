# Master Multiomics Simulation Benchmarking Pipeline

Orchestrates comprehensive benchmarking across all 7 evaluation
categories using strictly the user's provided inputs:

- `ref_multi`: List of reference omics matrices (e.g.
  `list(rna = ..., atac = ...)`).

- `sim_multi`: List of simulated omics matrices (e.g.
  `list(rna = ..., atac = ...)`).

- `cell_types`: Vector/factor of cell-type identities.

- `batch_info`: Vector/factor of batch or donor annotations.

- Computational resource usage: `cpu_time`, `memory_mb`, `system_time`.

- Trajectory features inferred directly from scRNA-seq counts.

## Usage

``` r
evaluate_multiomics_accuracy(
  ref_multi,
  sim_multi,
  cell_types = NULL,
  batch_info = NULL,
  cpu_time = NULL,
  memory_mb = NULL,
  system_time = NULL,
  elapsed_time = NULL,
  feature_pairs = NULL,
  compute_bivariate = FALSE,
  threads = 1,
  verbose = TRUE
)
```

## Arguments

- ref_multi:

  Named list of reference matrices for Modality 1 and 2.

- sim_multi:

  Named list of simulated matrices for Modality 1 and 2.

- cell_types:

  Optional factor or vector of cell type labels for cells.

- batch_info:

  Optional factor or vector of batch labels for cells.

- cpu_time:

  Optional numeric value of CPU execution time in seconds.

- memory_mb:

  Optional numeric value of peak memory allocation in MB.

- system_time:

  Optional numeric value of system CPU time in seconds.

- elapsed_time:

  Optional numeric value of wall-clock elapsed time in seconds.

- feature_pairs:

  Optional 2-column data.frame of linked feature pairs.

- compute_bivariate:

  Logical, whether to compute 2D bivariate tests. Default FALSE for
  speed.

- threads:

  Number of CPU threads. Default is 1.

- verbose:

  Logical, whether to print execution progress. Default is TRUE.

## Value

A list containing detailed results and a unified tidy master summary
table.
