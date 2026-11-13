# Compute dynverse Branch Assignment Overlap (F1_branches)

Evaluates recovery and relevance of cellular MST branch assignments
between reference and simulated trajectories via bidirectional Jaccard
overlap and harmonic mean F1 (matching
[`dyneval::calculate_mapping`](https://rdrr.io/pkg/dyneval/man/calculate_mapping.html)
`F1_branches`).

## Usage

``` r
calc_dynverse_f1_branches(ref_traj, sim_traj, nn_map = NULL)
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

Numeric `F1_branches` score in `[0, 1]` (higher is better).
