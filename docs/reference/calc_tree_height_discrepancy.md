# Calculate Lineage Tree Branch Height Discrepancy (RMSE)

Computes the Root Mean Squared Error (RMSE) between sorted branch
heights of hierarchical lineage trees. Accepts either precomputed
`hclust` objects OR scRNA-seq count matrices + `cell_types` (which will
be automatically inferred).

## Usage

``` r
calc_tree_height_discrepancy(
  ref_tree,
  sim_tree,
  cell_types_ref = NULL,
  cell_types_sim = NULL
)
```

## Arguments

- ref_tree:

  Hierarchical clustering (`hclust`) object from reference data, OR
  reference count matrix.

- sim_tree:

  Hierarchical clustering (`hclust`) object from simulated data, OR
  simulated count matrix.

- cell_types_ref:

  Optional cell type labels for reference (if count matrix provided).

- cell_types_sim:

  Optional cell type labels for simulation (if count matrix provided).

## Value

Root mean squared error between branch heights.
