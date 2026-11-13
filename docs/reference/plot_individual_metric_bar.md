# Plot Individual Metric Performance Barplot

Generates a publication-quality barplot comparing simulator methods on a
single evaluation metric, with exact score labels, direction awareness,
and performance ranking.

## Usage

``` r
plot_individual_metric_bar(
  benchmark_data,
  metric,
  score_type = c("normalized", "raw"),
  palette = NULL,
  base_size = 12
)
```

## Arguments

- benchmark_data:

  A benchmark data.frame or named list containing benchmark tables.

- metric:

  Character. The name of the metric to plot (e.g. `"KS"`,
  `"Wasserstein"`, `"ARI"`).

- score_type:

  Character. Either `"normalized"` (standardized \[0, 1\] fidelity
  score; default) or `"raw"` (original unnormalized metric value).

- palette:

  Optional named character vector of simulator colors.

- base_size:

  Numeric. Base font size. Default `12`.

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
