# Compute Weighted Feature Importance Correlation (wcor_features, Saelens et al., 2019)

Trains ranger Random Forest regression models to predict
cell-to-milestone geodesic distances from HVG expression values
(`num.trees = 10000`, `mtry_prop = 0.01` (1 percent of features),
`importance = "impurity"`), averages Gini impurity importance across
trajectory milestones for both reference and simulated trajectories, and
computes the reference-weighted Pearson correlation (`wcor_features`)
clipped to `[0, 1]`.

## Usage

``` r
calc_wcor_features(
  ref_traj,
  sim_traj,
  ref_expr,
  sim_expr,
  num_trees = 10000L,
  mtry_prop = 0.01,
  threads = 1L,
  seed = 42L
)
```

## Arguments

- ref_traj:

  Reference trajectory list containing `cell_milestone_geodesic`.

- sim_traj:

  Simulated trajectory list containing `cell_milestone_geodesic`.

- ref_expr:

  Reference log-normalized HVG matrix (`hvgs x cells_ref`).

- sim_expr:

  Simulated log-normalized HVG matrix (`hvgs x cells_sim`).

- num_trees:

  Number of Random Forest trees per milestone (default `10000L`, per
  Saelens et al. 2019).

- mtry_prop:

  Fraction of features sampled per split (default `0.01`, i.e., 1
  percent).

- threads:

  Number of threads for ranger (default `1L`).

- seed:

  Random seed for reproducible forest fitting (default `42L`).

## Value

A list containing `wcor_features` (weighted correlation clipped to
`[0, 1]`), `cor_features` (unweighted correlation clipped to `[0, 1]`),
`ref_importance`, and `sim_importance`.
