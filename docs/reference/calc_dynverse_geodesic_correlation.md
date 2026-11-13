# Compute dynverse Geodesic Distance Ordering Correlation

Computes pairwise cell-to-milestone geodesic distances along each
trajectory network, aligns cells across the shared HVG PCA manifold, and
computes the Spearman rank correlation clipped to `[0, 1]` (matching
`dyneval`'s `correlation` metric).

## Usage

``` r
calc_dynverse_geodesic_correlation(ref_traj, sim_traj, nn_map = NULL)
```

## Arguments

- ref_traj:

  Reference trajectory list.

- sim_traj:

  Simulated trajectory list.

- nn_map:

  Precomputed nearest-neighbor index map between reference and simulated
  cells.

## Value

Numeric geodesic distance correlation in `[0, 1]` (higher is better).
