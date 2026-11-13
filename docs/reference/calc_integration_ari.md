# Calculate Integration Adjusted Rand Index (ARI)

Computes the Adjusted Rand Index (ARI) evaluating the concordance
between clusters identified in the integrated/joint multiomics embedding
and ground-truth cell-type labels, following `scib`, `scCross`,
`scMultiSim`, `Matilda`, and `scMoMtF`.

## Usage

``` r
calc_integration_ari(embed, cell_type, n_clusters = NULL)
```

## Arguments

- embed:

  Numeric matrix of cell coordinates (cells x dimensions).

- cell_type:

  Vector or factor of ground-truth cell-type labels.

- n_clusters:

  Number of clusters (defaults to the number of unique cell types).

## Value

Numeric scalar between 0 and 1, where 1 indicates perfect clustering
match.

## References

Hubert & Arabie (1985) Comparing partitions. Journal of Classification.
Luecken et al. (2022) Benchmarking atlas-level data integration. Nature
Methods. Ding et al. (2023) scCross. Li et al. (2023) scMultiSim. Cell
Systems. Yang et al. (2023) Matilda. Nucleic Acids Research. Lan et al.
(2024) scMoMtF.

## Examples

``` r
set.seed(42)
emb <- matrix(rnorm(60 * 4), 60, 4)
ct <- rep(c("TypeA", "TypeB"), each = 30)
calc_integration_ari(emb, ct)
#> [1] 0
```
