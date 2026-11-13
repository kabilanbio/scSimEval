# Multi-Dimensional Scaling (MDS) Ordination of Evaluation Metrics or Simulators

Projects evaluation metric profiles or simulator performances into a 2D
MDS space for benchmarking single-cell simulation methods.

## Usage

``` r
plot_metric_mds(
  benchmark_data,
  ordination_by = c("simulators", "summaries"),
  by_category = FALSE,
  as_list = FALSE,
  category = NULL,
  exclude_categories = c("(VII) Trajectory & Lineage Topology",
    "(VIII) Computational Scalability & Efficiency"),
  metrics = NULL,
  normalize_scores = TRUE,
  palette = NULL,
  title = NULL,
  base_size = 11,
  ncol = 2
)
```

## Arguments

- benchmark_data:

  Benchmark summary data.frame, matrix, or list.

- ordination_by:

  Character. Either `"summaries"` (default; ordinated points represent
  biological and computational summaries) or `"simulators"` (ordinated
  points represent candidate simulators).

- by_category:

  Logical. If `TRUE` and `ordination_by = "summaries"`, generates a
  multi-panel grid displaying an MDS ordination for each evaluation
  category separately (skipping categories with fewer than 3 measures).
  Default is `FALSE`.

- as_list:

  Logical. If `TRUE` and `by_category = TRUE`, returns a named list of
  ggplot objects. Default `FALSE`.

- category:

  Optional character vector to filter by evaluation category (e.g.
  `"(I) Distributional Properties"`).

- exclude_categories:

  Character vector of categories to exclude from the ordination. Default
  is
  `c("(VII) Trajectory & Lineage Topology", "(VIII) Computational Scalability & Efficiency")`
  because these two categories contain only a few measures each, leaving
  the 6 canonical categories.

- metrics:

  Optional character vector of specific metrics to include.

- normalize_scores:

  Logical. If `TRUE` (default), scores are direction-aware normalized to
  \[0, 1\] so that metrics on disparate scales (e.g. MB vs distance) do
  not artificially distort Euclidean distance.

- palette:

  Optional named character vector of colors. For
  `ordination_by = "summaries"`, maps `"gene"`, `"cell"`, and
  `"global"`.

- title:

  Optional plot title.

- base_size:

  Numeric. Base font size. Default `11`.

- ncol:

  Integer. Number of columns when `by_category = TRUE`. Default `2`.

## Value

A `ggplot` or `patchwork` object.

## Details

When `ordination_by = "summaries"`, each point represents an evaluated
property/metric, colored by its biological level: `"gene"` (red),
`"cell"` (blue), or `"global"` (green). If `by_category = TRUE`, a
multi-panel dashboard is produced showing one MDS ordination per
category (by default across the 6 canonical categories with \>= 3
measures, excluding Trajectory and Scalability which each contain only 2
measures).

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
