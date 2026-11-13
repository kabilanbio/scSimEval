# Interactive Step-by-Step Zoom In and Zoom Out for Plots

Applies step-by-step magnification or custom coordinate window zooming
to any `ggplot2` figure produced by `scSimEval` (including comparative
bubble matrices, evaluation summaries, distributional QC, metric
heatmaps, and ordination plots).

## Usage

``` r
zoom_plot(p, zoom = 1, step = 0, xlim = NULL, ylim = NULL, center = NULL)
```

## Arguments

- p:

  A `ggplot` or `patchwork` object to zoom.

- zoom:

  Numeric zoom factor (\> 1.0 zooms in, \< 1.0 zooms out). Default is
  1.0 (no change).

- step:

  Integer step counter. Each positive step (+1) magnifies by +25%
  (1.25x), while each negative step (-1) demagnifies by -20% (0.8x).

- xlim:

  Optional numeric range `c(xmin, xmax)` to set exact horizontal limits.

- ylim:

  Optional numeric range `c(ymin, ymax)` to set exact vertical limits.

- center:

  Optional numeric vector `c(x_center, y_center)` specifying the focal
  point of zooming. Defaults to the center of the current plot
  coordinate span.

## Value

A modified `ggplot` or `patchwork` object with adjusted cartesian
coordinates.

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
p <- plot_evaluation_summary(res)
p_zoom_in <- zoom_plot(p, step = 1)
p_zoom_out <- zoom_plot(p, step = -1)
```
