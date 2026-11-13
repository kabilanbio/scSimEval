# Compute Simulator Method Performance Leaderboard

Computes the overall performance leaderboard ranking simulation methods
across all benchmark evaluation metrics. Metrics are first
direction-inverted and min-max standardized into fidelity scores in \[0,
1\] (where 1.0 represents best observed performance). Simulators are
then rank-ordered by mean overall fidelity score.

## Usage

``` r
compute_method_leaderboard(benchmark_data)
```

## Arguments

- benchmark_data:

  A benchmark data.frame (such as `demo$benchmark_summary_table` or
  output from
  [`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md))
  or a named list of benchmark result tables.

## Value

A `data.frame` with columns:

- `Overall_Rank`: Integer rank (1 = top-performing simulator).

- `Method`: Simulator method name.

- `Average_Fidelity`: Formatted percentage string (e.g., `"58.5%"`).

- `Fidelity_Score`: Numeric composite fidelity score in \[0, 1\] rounded
  to 4 decimals.

## Examples

``` r
data(example_scrna, package = "scSimEval")
res <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim)
#> [1/5] Extracting cell-level properties...
#> [2/5] Extracting feature-level properties...
#> [3/5] Computing univariate accuracy metrics...
#> [4/5] Computing bivariate metrics (Fasano-Franceschini, KDE zstat, 2D EMD)...
#> [5/5] Computing zero-probability and manifold distances...
#> Evaluating batch integration & confounder mixing...
#> Evaluating scRNA-seq biological signal fidelity (Libra Pseudobulk edgeR-QLF)...
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> [6/7] Evaluating trajectory and lineage topology (HVG + dynverse, scRNA-seq only)...
#> Unimodal accuracy evaluation complete.
lb <- compute_method_leaderboard(res)
```
