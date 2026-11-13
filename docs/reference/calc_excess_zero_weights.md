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

  Count matrix (features x cells).

## Value

A list with mean excess zero weight, gene-level weights, and estimated
zero-inflation rate.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_excess_zero_weights(example_scrna$ref)
#> $mean_excess_zero_weight
#> [1] 0.125175
#> 
#> $median_excess_zero_weight
#> [1] 0.05941929
#> 
#> $mean_zero_inflation
#> [1] 0.02386441
#> 
#> $gene_excess_weights
#>  [1] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.09777775
#>  [7] 0.40423218 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> [13] 0.00000000 0.44714706 0.08562163 0.18806978 0.27267149 0.06694220
#> [19] 0.00000000 0.00000000 0.01533320 0.00000000 0.52388338 0.00000000
#> [25] 0.00000000 0.00000000 0.00000000 0.00000000 0.49862464 0.28463070
#> [31] 0.00000000 0.02307153 0.31073986 0.00000000 0.06367513 0.15088815
#> [37] 0.24689382 0.00000000 0.09061268 0.31942074 0.21672462 0.18998773
#> [43] 0.07329427 0.00000000 0.22786326 0.00000000 0.34310829 0.06448393
#> [49] 0.05516346 0.00000000 0.30971601 0.19499308 0.28160864 0.37442912
#> [55] 0.00000000 0.27711280 0.26124999 0.45666487 0.09386449 0.00000000
#> 
```
