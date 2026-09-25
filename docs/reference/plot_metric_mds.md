# Multi-Dimensional Scaling (MDS) Ordination of Evaluation Metrics or Simulators

Projects evaluation metric profiles or simulator performances into a 2D
MDS space for benchmarking single-cell simulation methods.

## Usage

``` r
plot_metric_mds(
  benchmark_data,
  ordination_by = c("summaries", "simulators"),
  by_category = FALSE,
  as_list = FALSE,
  category = NULL,
  exclude_categories = c("(VI) Trajectory & Lineage Dynamics",
    "(VIII) Computational Scalability"),
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

- category:

  Optional character vector to filter by evaluation category (e.g.
  `"(I) Distributional Properties"`).

- exclude_categories:

  Character vector of categories to exclude from the ordination. Default
  is
  `c("(VI) Trajectory & Lineage Dynamics", "(VIII) Computational Scalability")`
  because these two categories contain only two measures each, leaving
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
