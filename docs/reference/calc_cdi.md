# Compute Clustering Deviation Index (CDI) for Single-Cell Clustering Evaluation

Evaluates the unsupervised goodness-of-fit of single-cell clustering
labels directly on count matrices using the official Bioconductor CDI
package (Fang et al., *Genome Biology*, 2022; DOI:
10.1186/s13059-022-02825-5). CDI models UMI counts under cluster- and
gene-specific Negative Binomial distributions adjusted by cell size
factors, returning AIC- and BIC-penalized deviation scores.

## Usage

``` r
calc_cdi(
  counts,
  cluster_labels,
  batch_label = NULL,
  size_factors = NULL,
  nfeature = 500L,
  method = c("wds", "vst"),
  zp_threshold = 0.95,
  lrt_pval_threshold = 0.01,
  max_cells = 1000L,
  normalize_per_obs = TRUE,
  threads = 1L,
  seed = 42L,
  top_features = NULL
)
```

## Arguments

- counts:

  Count matrix (`matrix`, sparse `dgCMatrix`, or `SingleCellExperiment`;
  genes x cells or cells x genes).

- cluster_labels:

  Vector or factor of cluster / cell-type assignments for each cell.

- batch_label:

  Optional vector or factor of batch assignments for each cell (enables
  CDI's batch-common vs. batch-specific NB likelihood ratio test).

- size_factors:

  Optional numeric vector of positive cell-specific size factors. If
  `NULL`, computed via
  [`CDI::size_factor()`](https://rdrr.io/pkg/CDI/man/size_factor.html).

- nfeature:

  Integer number of informative feature genes to select via
  [`CDI::feature_gene_selection()`](https://rdrr.io/pkg/CDI/man/feature_gene_selection.html)
  (default `500L`, per Fang et al., 2022).

- method:

  Feature selection method passed to
  [`CDI::feature_gene_selection()`](https://rdrr.io/pkg/CDI/man/feature_gene_selection.html):
  `"wds"` (Working Dispersion Score, default) or `"vst"`.

- zp_threshold:

  Numeric zero-proportion cutoff in `(0, 1)` used to pre-filter
  ultra-sparse features prior to WDS selection (default `0.95`).
  Automatically relaxed for ultra-sparse modalities such as scATAC-seq.

- lrt_pval_threshold:

  Numeric p-value threshold in `(0, 1)` for multi-batch LRT in
  [`CDI::calculate_CDI()`](https://rdrr.io/pkg/CDI/man/calculate_CDI.html)
  (default `0.01`).

- max_cells:

  Maximum number of cells evaluated via proportional stratified
  subsampling across clusters (default `1000L`; set `Inf` or `NULL` to
  evaluate all cells).

- normalize_per_obs:

  Logical; if `TRUE` (default), divides raw `CDI_AIC` and `CDI_BIC` by
  \\N\_{\mathrm{eval}} \times G\_{\mathrm{selected}}\\ so scores are
  scale-invariant across simulated datasets of differing dimensions.

- threads:

  Integer number of parallel worker threads for `BiocParallel` (default
  `1L`).

- seed:

  Integer random seed used only during stratified cell subsampling
  (default `42L`).

- top_features:

  Optional deprecated alias for `nfeature`.

## Value

A list containing:

- cdi_aic:

  Per-observation normalized (or raw if `normalize_per_obs = FALSE`)
  CDI-AIC score (lower indicates better clustering, sensitive to
  subtypes).

- cdi_bic:

  Per-observation normalized (or raw if `normalize_per_obs = FALSE`)
  CDI-BIC score (lower indicates better clustering, favors main cell
  types).

- cdi_aic_raw:

  Unnormalized `CDI_AIC` returned directly by
  [`CDI::calculate_CDI()`](https://rdrr.io/pkg/CDI/man/calculate_CDI.html).

- cdi_bic_raw:

  Unnormalized `CDI_BIC` returned directly by
  [`CDI::calculate_CDI()`](https://rdrr.io/pkg/CDI/man/calculate_CDI.html).

- loglik:

  Total maximized NB log-likelihood (`-neg_llk_val`) across selected
  features.

- n_clusters:

  Number of clusters evaluated.

- n_features:

  Number of selected feature genes evaluated.

- n_cells:

  Number of cells evaluated.

## Details

To minimize time and space complexity from \\O(G \cdot N)\\ dense
allocations down to \\O(\mathrm{nnz}(X) + G_0 \cdot
N\_{\mathrm{eval}})\\, sparse `dgCMatrix` inputs are
stratified-subsampled across cells and pre-filtered by zero-proportion
in sparse space prior to Working Dispersion Score (WDS) feature
selection
([`CDI::feature_gene_selection`](https://rdrr.io/pkg/CDI/man/feature_gene_selection.html)),
median-of-ratios size factor estimation
([`CDI::size_factor`](https://rdrr.io/pkg/CDI/man/size_factor.html)),
and Negative Binomial MLE evaluation
([`CDI::calculate_CDI`](https://rdrr.io/pkg/CDI/man/calculate_CDI.html)).

## Examples

``` r
counts <- matrix(stats::rpois(1200, lambda = 3), nrow = 30, ncol = 40)
cl <- factor(rep(c("TypeA", "TypeB"), each = 20))
calc_cdi(counts, cl)
#> $cdi_aic
#> [1] 4.053388
#> 
#> $cdi_bic
#> [1] 4.222276
#> 
#> $cdi_aic_raw
#> [1] 4864.066
#> 
#> $cdi_bic_raw
#> [1] 5066.731
#> 
#> $loglik
#> [1] -2312.033
#> 
#> $n_clusters
#> [1] 2
#> 
#> $n_features
#> [1] 30
#> 
#> $n_cells
#> [1] 40
#> 
```
