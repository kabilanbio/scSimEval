# Plot Cross-Modal Regulatory Coupling (RNA \<-\> ATAC)

For multiomics benchmarking: compares Spearman cross-modality
correlations between linked gene expression (scRNA-seq) and chromatin
accessibility (scATAC-seq) pairs.

## Usage

``` r
plot_cross_modal_coupling(
  ref_rna,
  ref_atac,
  sim_rna,
  sim_atac,
  n_features = 20,
  palette = c("#2E86AB", "#E74C3C")
)
```

## Arguments

- ref_rna:

  Reference scRNA-seq matrix (genes x cells).

- ref_atac:

  Reference scATAC-seq matrix (peaks x cells).

- sim_rna:

  Simulated scRNA-seq matrix (genes x cells).

- sim_atac:

  Simulated scATAC-seq matrix (peaks x cells).

- n_features:

  Integer. Number of linked gene-peak pairs to evaluate. Default `20`.

- palette:

  Character vector of 2 colors. Default `c("#2E86AB", "#E74C3C")`.

## Value

A `ggplot` object of cross-modal correlation bars.

## Examples

``` r
data(example_multiomics)
p <- plot_cross_modal_coupling(
  example_multiomics$ref_multi$rna, example_multiomics$ref_multi$atac,
  example_multiomics$sim_multi$rna, example_multiomics$sim_multi$atac
)
if (requireNamespace("ggplot2", quietly = TRUE)) print(p)
```
