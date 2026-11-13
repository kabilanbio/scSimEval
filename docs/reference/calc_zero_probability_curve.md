# Fit Zero-Probability Dropout Curve (Splatter & ZINB-WaVE)

Fits an empirical logistic dropout curve: logit(P(Y=0)) = beta_0 +
beta_1 \* log(mu) to model the scRNA-seq dropout relationship with mean
expression.

## Usage

``` r
calc_zero_probability_curve(counts)
```

## Arguments

- counts:

  Count matrix (genes x cells) or SingleCellExperiment.

## Value

A list containing intercept, slope, midpoint (inflection point), and
R-squared.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_zero_probability_curve(example_scrna$ref)
#> $intercept
#> [1] 0.1664043
#> 
#> $slope
#> [1] -1.259067
#> 
#> $midpoint
#> [1] 0.1321648
#> 
#> $r_squared
#> [1] 0.7142614
#> 
```
