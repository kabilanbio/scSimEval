# Calculate Homogeneity, Completeness, and V-Measure

Integrated from scCluBench (Xu et al., AAAI 2026).

## Usage

``` r
calc_homogeneity_completeness_v_measure(pred, truth, beta = 1)
```

## Arguments

- pred:

  Vector of predicted cluster labels.

- truth:

  Vector of ground truth cell type labels.

- beta:

  Weight of completeness vs. homogeneity (default 1).

## Value

A named list containing homogeneity, completeness, and v_measure.
