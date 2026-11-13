# Run edgeR QLF Differential Expression with Cellular Detection Rate (edgeRQLFDetRate)

Implements the benchmark-optimal `edgeRQLFDetRate` algorithm (Soneson &
Robinson, 2018; Duo et al., 2024 / Simpipe) for single-cell differential
expression analysis. Fits a quasi-likelihood negative binomial
generalized linear model including the cellular detection rate (DetRate:
fraction of expressed genes per cell) as a covariate.

## Usage

``` r
run_edger_qlf_detrate(
  counts,
  group,
  cdr = NULL,
  fdr_cutoff = 0.05,
  logfc_cutoff = 0,
  min_cells = 3
)
```

## Arguments

- counts:

  Single-cell count matrix (genes x cells), or SingleCellExperiment /
  Seurat object.

- group:

  Factor or vector specifying the comparison cell groups (must contain
  exactly 2 levels for pairwise DE).

- cdr:

  Optional numeric vector of cellular detection rates per cell. If NULL,
  computed automatically as fraction of non-zero genes.

- fdr_cutoff:

  FDR / adjusted p-value threshold for significance (default 0.05).

- logfc_cutoff:

  Absolute log2 fold-change threshold (default 0).

- min_cells:

  Minimum number of cells expressing a gene to test (default 3).

## Value

A data.frame with columns:

- gene:

  Gene identifier

- logFC:

  log2 fold change

- logCPM:

  log2 counts per million

- F_stat:

  Quasi-likelihood F-statistic (or t-statistic in fallback)

- PValue:

  Raw p-value

- FDR:

  Benjamini-Hochberg adjusted p-value

- is_deg:

  Logical flag indicating whether the gene meets the significance
  threshold

## Examples

``` r
data(example_scrna, package = "scSimEval")
top2 <- names(sort(table(example_scrna$cell_types), decreasing = TRUE))[1:2]
idx <- which(example_scrna$cell_types %in% top2)
res <- run_edger_qlf_detrate(example_scrna$ref[, idx], example_scrna$cell_types[idx])
#> calcNormFactors has been renamed to normLibSizes
head(res)
#>            gene       logFC   logCPM     F_stat     PValue       FDR is_deg
#> Gene_01 Gene_01  0.06639992 15.04487 0.04541237 0.83125697 0.9236189  FALSE
#> Gene_02 Gene_02  0.63252452 14.58323 3.77514999 0.05207991 0.3124794  FALSE
#> Gene_03 Gene_03  0.69780444 14.79379 4.75012284 0.02934654 0.3124794  FALSE
#> Gene_04 Gene_04  0.30634068 14.69534 0.89011764 0.34549476 0.6909895  FALSE
#> Gene_05 Gene_05  0.54739063 15.07996 2.95097305 0.08589288 0.3594412  FALSE
#> Gene_06 Gene_06 -0.09283335 14.65058 0.08278491 0.77357118 0.8925821  FALSE
```
