# Plot Trajectory Dynamics Comparison (Pseudotime Q-Q Alignment)

Compares continuous differentiation pseudotime distributions between
reference and simulated scRNA-seq data via quantile-quantile (Q-Q)
alignment and overlaid KDE curves. Pseudotime is auto-inferred from
count matrices if not provided as precomputed vectors.

## Usage

``` r
plot_trajectory_comparison(
  ref_data,
  sim_data,
  palette = c("#2E86AB", "#E74C3C"),
  n_quantile_pts = 200
)
```

## Arguments

- ref_data:

  Reference count matrix (genes x cells) or precomputed pseudotime
  vector.

- sim_data:

  Simulated count matrix (genes x cells) or precomputed pseudotime
  vector.

- palette:

  Character vector of 2 colors. Default `c("#2E86AB", "#E74C3C")`.

- n_quantile_pts:

  Integer. Quantile points for Q-Q plot. Default `200`.

## Value

A `ggplot` object (or patchwork object if patchwork is available).

## Examples

``` r
data(example_scrna)
p <- plot_trajectory_comparison(example_scrna$ref, example_scrna$sim)
if (requireNamespace("ggplot2", quietly = TRUE)) print(p)
#> `geom_smooth()` using formula = 'y ~ x'
```
