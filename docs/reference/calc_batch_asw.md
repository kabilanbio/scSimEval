# Calculate Batch Average Silhouette Width (Batch ASW)

Evaluates batch mixing within each cell type in integrated or joint
multiomics embedding space, following the definitions in `scib` and
`scMultiSim`. For each cell type, the silhouette width of cells with
respect to batch labels is computed (\\-1 \le s_i \le 1\\), where values
near 0 indicate well-mixed batches. The per-cell-type score is \\1 -
\|s_i\|\\, averaged across all cells and cell types.

## Usage

``` r
calc_batch_asw(embed, batch, cell_type = NULL)
```

## Arguments

- embed:

  Numeric matrix of cell coordinates (cells x dimensions, e.g. PCA or
  joint embedding).

- batch:

  Vector or factor of batch or modality labels for each cell.

- cell_type:

  Optional vector or factor of cell-type labels for each cell.

## Value

Numeric scalar between 0 and 1, where 1 indicates optimal batch mixing.

## References

Luecken et al. (2022) Benchmarking atlas-level data integration in
single-cell genomics. Nature Methods. Li et al. (2023) scMultiSim:
simulation of multi-modality single cell data guided by gene regulatory
networks. Cell Systems.

## Examples

``` r
set.seed(42)
emb <- matrix(rnorm(100 * 10), 100, 10)
bt <- rep(c("B1", "B2"), each = 50)
ct <- rep(c("T1", "T2"), 50)
calc_batch_asw(emb, bt, ct)
#> [1] 0.9562664
```
