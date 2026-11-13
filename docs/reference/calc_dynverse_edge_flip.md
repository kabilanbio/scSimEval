# Compute dynverse Edge Flip Trajectory Topology Similarity

Quantifies topological similarity between reference and simulated
milestone graphs via minimum edge flip edit distance
([`dyneval::calculate_edge_flip`](https://rdrr.io/pkg/dyneval/man/calculate_edge_flip.html)).

## Usage

``` r
calc_dynverse_edge_flip(ref_traj, sim_traj)
```

## Arguments

- ref_traj:

  Reference trajectory list returned by
  [`infer_dynverse_trajectory`](https://kabilanbio.github.io/scSimEval/reference/infer_dynverse_trajectory.md).

- sim_traj:

  Simulated trajectory list returned by
  [`infer_dynverse_trajectory`](https://kabilanbio.github.io/scSimEval/reference/infer_dynverse_trajectory.md).

## Value

Numeric edge-flip similarity in `[0, 1]` (higher is better).
