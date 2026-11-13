# Plot Benchmark Metric Performance Barplots by Category

Generates publication-ready faceted barplots comparing simulator methods
across all individual evaluation metrics belonging to a specified
benchmark category. Each metric is presented in its own sub-panel with
exact numeric score labels, direction-awareness indicators (indicating
whether higher or lower values are optimal), and standardized or
original raw score formulation.

## Usage

``` r
plot_category_metric_bars(
  benchmark_data,
  category,
  score_type = c("normalized", "raw"),
  palette = NULL,
  ncol = NULL,
  base_size = 11
)
```

## Arguments

- benchmark_data:

  A benchmark data.frame (e.g. `demo$benchmark_summary_table` or output
  from
  [`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md))
  or a named list of benchmark result tables.

- category:

  Character. The benchmark category to visualize (e.g.,
  `"(I) Distributional Properties"` or `"Distributional Properties"`).

- score_type:

  Character. Either `"normalized"` (standardized \[0, 1\] fidelity score
  where 1.0 is optimal; default) or `"raw"` (unnormalized original
  metric values in native measurement units).

- palette:

  Optional named character vector of simulator colors.

- ncol:

  Integer. Optional number of columns in facet layout. Defaults to
  dynamic layout based on metric count.

- base_size:

  Numeric. Base font size. Default `11`.

## Value

A `ggplot` object.

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
```
