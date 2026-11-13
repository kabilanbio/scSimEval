# Plot Single-Cell Summary Distribution Quality (Comparative Overlays)

Compares empirical reference and simulated datasets across 18
distributional properties: 7 univariate cell-level properties, 7
univariate feature-level properties, and 4 bivariate property pairs,
either as empirical & joint density overlays or as quantile-quantile
(Q-Q) concordance plots.

## Usage

``` r
plot_distribution_qc(
  ref_data,
  sim_data,
  layout = c("density", "qq", "comprehensive"),
  properties = NULL,
  title = "Data properties",
  ref_name = "Original",
  sim_name = "Simulation",
  palette = NULL,
  cell_types = NULL,
  base_size = 9
)
```

## Arguments

- ref_data:

  Numeric count matrix for empirical reference (features x cells).

- sim_data:

  Numeric count matrix for simulated data (features x cells), or a named
  list of simulated count matrices representing multiple simulators
  (e.g. `list("scDesign3" = m1, "Splatter" = m2)`).

- layout:

  Character. Visualization layout: `"qq"` (18-panel quantile-quantile
  diagnostic plots comparing real vs. simulated quantiles on a common
  `n = 1024` grid with KS, MAD, MAE, RMSE, and PCC annotations) or
  `"density"` / `"comprehensive"` (18-panel empirical univariate density
  curves and bivariate joint distributions across the exact same 18
  properties).

- properties:

  Optional character vector of specific properties or group keywords
  (`"all"`, `"cell"`, `"feature"`, `"bivariate"`, or a single property
  key) to include.

- title:

  Character. Main plot title. Default `"Data properties"` (uses
  `"Real vs. Simulated Q-Q Plots"` for `layout = "qq"` and
  `"Real vs. Simulated Empirical & Joint Density Plots"` for density
  layout when left at default). Set to `NULL` for no title.

- ref_name:

  Character. Label for the empirical reference in the legend. Default
  `"Original"`.

- sim_name:

  Character. Label for the simulation if a single matrix is provided.
  Default `"Simulation"`.

- palette:

  Optional named character vector of colors. Default uses warm red for
  reference and classic steel blue for simulation.

- cell_types:

  Optional factor or character vector (retained for backward
  compatibility).

- base_size:

  Numeric. Base font size. Default `9`.

## Value

A `ggplot` or `patchwork` composite object containing the comparative
panels.

## Details

![Comparative Single-Cell Distribution
QC](figures/comparative_distribution_qc.png)

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
data(example_scrna)
p <- plot_distribution_qc(example_scrna$ref, example_scrna$sim)
if (requireNamespace("ggplot2", quietly = TRUE)) print(p)
```
