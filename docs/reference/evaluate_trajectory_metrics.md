# Full Trajectory Accuracy Evaluation

Evaluates trajectory preservation using pseudotime correlation, branch
height discrepancy, and univariate distribution distance metrics on
pseudotime. If raw count matrices are provided, pseudotime and lineage
trees are automatically inferred directly from the scRNA-seq data.

## Usage

``` r
evaluate_trajectory_metrics(
  ref_data,
  sim_data,
  cell_types_ref = NULL,
  cell_types_sim = NULL
)
```

## Arguments

- ref_data:

  Reference count matrix (genes x cells) OR numeric vector of
  pseudotimes.

- sim_data:

  Simulated count matrix (genes x cells) OR numeric vector of
  pseudotimes.

- cell_types_ref:

  Optional vector of reference cell types (for lineage tree inference).

- cell_types_sim:

  Optional vector of simulated cell types (for lineage tree inference).

## Value

A named list of trajectory accuracy metrics.
