# Calculate Proportion of Outliers in a Numeric Vector

Uses the standard Tukey 1.5 \* IQR criterion on quantiles, as used in
simpipe.

## Usage

``` r
calc_outlier_proportion(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Outlier proportion (between 0 and 1).

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_outlier_proportion(ref)
#> [1] 0
```
