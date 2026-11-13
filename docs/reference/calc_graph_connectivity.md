# Calculate KNN Graph Connectivity

Evaluates whether cells of the same cell type form a connected component
in the k-nearest neighbor (k-NN) graph across modalities/batches,
following `scCross`, `scMultiSim`, and `scib`. For each cell type, the
relative size of the largest connected component (LCC) in the induced
subgraph is calculated: \$\|LCC\| / \|V_c\|\$. The overall score is the
mean across all cell types.

## Usage

``` r
calc_graph_connectivity(embed, cell_type, k = 15)
```

## Arguments

- embed:

  Numeric matrix of cell coordinates (cells x dimensions).

- cell_type:

  Vector or factor of cell-type labels for each cell.

- k:

  Number of nearest neighbors to construct the graph (default 15).

## Value

Numeric scalar between 0 and 1, where 1 indicates complete graph
connectivity.

## References

Ding et al. (2023) scCross: deep generative modeling for single-cell
cross-modality integration. Li et al. (2023) scMultiSim: simulation of
multi-modality single cell data. Cell Systems. Luecken et al. (2022)
Benchmarking atlas-level data integration. Nature Methods.

## Examples

``` r
set.seed(42)
emb <- matrix(rnorm(60 * 5), 60, 5)
ct <- rep(c("A", "B"), each = 30)
calc_graph_connectivity(emb, ct, k = 10)
#> [1] 1
```
