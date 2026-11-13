# Calculate Trajectory Graph Alignment Score

Evaluates trajectory topology preservation by measuring the concordance
of pairwise graph distances between landmark cells along the inferred
trajectory.

## Usage

``` r
calc_trajectory_graph_alignment(
  ref_data,
  sim_data,
  ref_pseudotime = NULL,
  sim_pseudotime = NULL,
  n_landmarks = 50
)
```

## Arguments

- ref_data:

  Reference expression/count matrix OR pseudotime vector.

- sim_data:

  Simulated expression/count matrix OR pseudotime vector.

- ref_pseudotime:

  Optional precomputed reference pseudotime.

- sim_pseudotime:

  Optional precomputed simulated pseudotime.

- n_landmarks:

  Number of trajectory landmark points (default 50).

## Value

Numeric trajectory graph alignment score in \[0, 1\] (higher is better).

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_trajectory_graph_alignment(example_scrna$ref, example_scrna$sim)
#> [1] 0.8718393
```
