# Principal Component Analysis (PCA) Dashboard of Benchmark Metrics and Methods

Generates a dual-panel PCA dashboard matching single-cell simulation
benchmarking publication standards (e.g. Crowell et al., Nature
Biotechnology; Soneson et al., Genome Biology).

## Usage

``` r
plot_metric_pca(
  benchmark_data,
  panel = c("both", "methods", "loadings"),
  by_category = FALSE,
  as_list = FALSE,
  category = NULL,
  exclude_categories = NULL,
  metrics = NULL,
  palette = NULL,
  top_n_loadings = 14,
  base_size = 11
)
```

## Arguments

- benchmark_data:

  Benchmark summary data.frame, matrix, or list.

- panel:

  Character. Which panel(s) to render: `"both"` (default; 2-panel
  stacked layout), `"methods"` (panel a only), or `"loadings"` (panel b
  only).

- category:

  Optional character vector to filter by evaluation category.

- exclude_categories:

  Optional character vector of categories to exclude (e.g. scalability).

- metrics:

  Optional character vector of specific metrics to include.

- palette:

  Optional named character vector of colors for simulators.

- top_n_loadings:

  Integer or `NULL`. Maximum number of loading vectors to plot in panel
  b (ordered by vector magnitude in PC1-PC2 space). Default is `14` to
  prevent visual clutter; set to `NULL` or `62` to display all evaluated
  metrics.

- base_size:

  Numeric. Base font size. Default `11`.

## Value

A `ggplot` or `patchwork` object representing the ordination figure.

## Details

All metrics present in `benchmark_data` (e.g. all 62 canonical metrics
across categories) are mathematically considered in the PCA
decomposition. Panel `"both"` displays:

- **Panel a (Methods Ordination)**: Simulators mapped into PC1 vs PC2
  coordinates with solid black crosshairs, centroid markers, and white
  rounded badge labels with colored borders.

- **Panel b (PC Loadings Vectors)**: Directional loading vectors
  (arrows) radiating from the origin `(0, 0)` to loading coordinates
  `(w1, w2)`, colored by summary type (`"gene"` = red, `"cell"` = blue,
  `"global"` = green) with rounded badge labels. Pass
  `top_n_loadings = NULL` or `62` to render all metric loading vectors.
