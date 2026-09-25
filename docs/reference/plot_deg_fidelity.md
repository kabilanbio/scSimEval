# Plot Differentially Expressed Gene (DEG) Fidelity

Produces a two-panel 600 DPI visualization evaluating differential
expression fidelity between empirical reference and simulated datasets:

- **Panel A (Effect Size Concordance):** Scatter plot comparing
  reference vs simulated log2 fold-changes across all features,
  annotated with Pearson correlation (\$r\$) and colored by DEG
  discovery status (Shared DEG, Reference Only, Simulation Only, Not
  DE).

- **Panel B (Multi-Framework Metric Benchmark):** Horizontal bar chart
  displaying fidelity metrics spanning Simpipe, SimBench, and Shaky
  Foundations.

## Usage

``` r
plot_deg_fidelity(deg_res, base_size = 11)
```

## Arguments

- deg_res:

  Result object returned by
  [`evaluate_deg_fidelity`](https://kabilanbio.github.io/scSimEval/reference/evaluate_deg_fidelity.md).

- base_size:

  Base font size for ggplot2 elements (default 11).

## Value

A ggplot object.
