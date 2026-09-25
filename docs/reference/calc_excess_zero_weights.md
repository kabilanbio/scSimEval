# Calculate Posterior Excess Zero Weights (zingeR & ZINB-WaVE)

Calculates the posterior probability of zero counts being technical
dropouts (excess zeros) versus biological sampling zeros under a
Negative Binomial model.

## Usage

``` r
calc_excess_zero_weights(counts)
```

## Arguments

- counts:

  Count matrix (genes x cells).

## Value

A list with mean excess zero weight, gene-level weights, and estimated
zero-inflation rate.
