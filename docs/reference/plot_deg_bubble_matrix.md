# Plot Multi-Framework Differential Expression Bubble Matrix

Renders a 600 DPI bubble matrix focused specifically on comparing
multiple single-cell simulation methods across the 15 differential
expression and biological signal metrics from Simpipe, SimBench, and
Shaky Foundations.

## Usage

``` r
plot_deg_bubble_matrix(
  data,
  methods_order = NULL,
  metrics_order = NULL,
  bubble_size_range = c(2, 8.5),
  base_size = 10,
  title = "Multi-Framework Differential Expression & Biological Signal Fidelity",
  subtitle =
    "Benchmarking simulation methods across Simpipe, SimBench, and Shaky Foundations"
)
```

## Arguments

- data:

  A list of result objects from
  [`evaluate_deg_fidelity`](https://kabilanbio.github.io/scSimEval/reference/evaluate_deg_fidelity.md)
  (e.g., `list("scDesign3" = res1, "Splatter" = res2)`) or a
  consolidated `data.frame` containing columns `Method`, `Framework`,
  `Metric`, and `Value`.

- methods_order:

  Optional character vector specifying custom method order.

- metrics_order:

  Optional character vector specifying custom metric order.

- bubble_size_range:

  Range of bubble radii (default `c(2, 8.5)`).

- base_size:

  Base font size (default 10).

- title:

  Plot title.

- subtitle:

  Plot subtitle.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming deg_res1 and deg_res2 are results from evaluate_deg_fidelity()
p_bubble <- plot_deg_bubble_matrix(
  data = list("scDesign3" = deg_res1, "Splatter" = deg_res2)
)
print(p_bubble)
} # }
```
