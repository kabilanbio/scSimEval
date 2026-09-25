# Calculate Fowlkes-Mallows Index (FMI)

Integrated from scCluBench (Xu et al., AAAI 2026). Computes the
geometric mean of pairwise cluster precision and recall.

## Usage

``` r
calc_fmi(pred, truth)
```

## Arguments

- pred:

  Vector of predicted cluster labels.

- truth:

  Vector of ground truth cell type labels.

## Value

FMI score between 0 and 1.
