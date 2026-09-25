# Plot Single-Cell Simulator Evaluation Summary

Produces a 600 DPI horizontal bar matrix ranking single-cell simulators
across the 8 canonical evaluation categories and overall composite
performance, styled in accordance with premier benchmark literature
(e.g., Nature Methods / Cell).

## Usage

``` r
plot_evaluation_summary(
  data,
  method_classes = NULL,
  category_colors = NULL,
  normalize_scores = TRUE,
  simulators_order = NULL,
  show_labels = TRUE,
  title = "Single-Cell Simulator Evaluation Summary",
  subtitle =
    "Average fidelity scores across 8 canonical evaluation categories and overall benchmark performance",
  base_size = 11
)

plot_benchmark_summary_bars(
  data,
  method_classes = NULL,
  category_colors = NULL,
  normalize_scores = TRUE,
  simulators_order = NULL,
  show_labels = TRUE,
  title = "Single-Cell Simulator Evaluation Summary",
  subtitle =
    "Average fidelity scores across 8 canonical evaluation categories and overall benchmark performance",
  base_size = 11
)

plot_summary_bars(
  data,
  method_classes = NULL,
  category_colors = NULL,
  normalize_scores = TRUE,
  simulators_order = NULL,
  show_labels = TRUE,
  title = "Single-Cell Simulator Evaluation Summary",
  subtitle =
    "Average fidelity scores across 8 canonical evaluation categories and overall benchmark performance",
  base_size = 11
)
```

## Arguments

- data:

  A benchmark summary data frame, a list containing
  `benchmark_summary_table`, or the output of `run_benchmark_suite`.

- method_classes:

  Optional named vector mapping method names to class/category labels
  (e.g., `c("Splatter" = "Class 1", "scDesign3" = "Class 2")`). If
  provided, rows are grouped into class panels.

- category_colors:

  Optional named vector specifying custom colors for category panels and
  overall score.

- normalize_scores:

  Logical indicating whether to normalize scores to \[0, 1\]
  (direction-aware where larger is better). Default is `TRUE`.

- simulators_order:

  Optional vector specifying a custom ordering of simulators. Default is
  ranking by descending Overall Score.

- show_labels:

  Logical indicating whether to display numerical score labels on/beside
  each bar. Default is `TRUE`.

- title:

  Character plot title. Default is
  `"Single-Cell Simulator Evaluation Summary"`.

- subtitle:

  Character plot subtitle.

- base_size:

  Numeric base font size (default 11).

## Value

A `ggplot` object representing the multi-panel evaluation summary.

## Details

![Single-Cell Simulator Evaluation
Summary](figures/evaluation_summary_bars.png)

## Examples

``` r
if (FALSE) { # \dontrun{
demo <- readRDS(system.file("shiny/scSimEvalApp/data/demo_benchmark_data.rds", package = "scSimEval"))
p <- plot_evaluation_summary(demo$benchmark_summary_table)
print(p)
} # }
```
