# Infer Docker-Free dynwrap Milestone Trajectory from Normalized Expression

Constructs a continuous Minimum Spanning Tree (MST) milestone trajectory
directly from log-normalized expression counts (typically HVG-filtered
via
[`preprocess_trajectory_hvgs`](https://kabilanbio.github.io/scSimEval/reference/preprocess_trajectory_hvgs.md)),
projects cells continuously onto MST edges, and computes pairwise
cell-to-milestone geodesic distances along the milestone network.

## Usage

``` r
infer_dynverse_trajectory(
  expr_mat,
  cell_types = NULL,
  n_pcs = 15L,
  pca_coords = NULL,
  milestone_prefix = "M"
)
```

## Arguments

- expr_mat:

  Log-normalized expression matrix (`genes x cells`).

- cell_types:

  Optional character vector or factor of cell-type identities.

- n_pcs:

  Number of PCA components for MST construction (default `15L`).

- pca_coords:

  Optional precomputed cell PCA matrix (`cells x dims`).

- milestone_prefix:

  Prefix for unsupervised cluster milestones when `cell_types` is
  `NULL`.

## Value

A list (augmented `dynwrap` trajectory object when dynwrap is available)
containing `milestone_ids`, `milestone_network`, `progressions`,
`milestone_percentages`, `pca_coords`, and `cell_milestone_geodesic`.

## Examples

``` r
data(example_scrna, package = "scSimEval")
prep <- preprocess_trajectory_hvgs(example_scrna$ref, cell_types_ref = example_scrna$cell_types)
traj <- infer_dynverse_trajectory(prep$ref_expr, prep$cell_types_ref)
#> Error in dimnames(x) <- dn: length of 'dimnames' [1] not equal to array extent
```
