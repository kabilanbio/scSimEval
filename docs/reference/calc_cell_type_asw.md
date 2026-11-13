# Calculate Cell-Type Average Silhouette Width (Cell-Type ASW)

Evaluates biological conservation and cell-type cluster separation in
integrated or joint multiomics space, following `scCross` and `scib`.
The silhouette width of cells with respect to ground-truth cell-type
labels is computed and scaled to \$\[0, 1\]\$ via \$(ASW + 1) / 2\$.

## Usage

``` r
calc_cell_type_asw(embed, cell_type)
```

## Arguments

- embed:

  Numeric matrix of cell coordinates (cells x dimensions).

- cell_type:

  Vector or factor of cell-type labels for each cell.

## Value

Numeric scalar between 0 and 1, where 1 indicates distinct, tight
cell-type clusters.

## References

Ding et al. (2023) scCross: deep generative modeling for single-cell
cross-modality integration and perturbation analysis. Luecken et al.
(2022) Benchmarking atlas-level data integration in single-cell
genomics. Nature Methods.

## Examples

``` r
set.seed(42)
emb <- matrix(rnorm(100 * 10), 100, 10)
ct <- rep(c("T1", "T2"), each = 50)
calc_cell_type_asw(emb, ct)
#> [1] 0.5009777
```
