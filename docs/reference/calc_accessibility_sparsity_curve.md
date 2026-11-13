# Fit Chromatin Accessibility-Sparsity Polynomial Curve (simATAC)

Fits a polynomial curve relating peak/bin mean accessibility to non-zero
cell proportion (NZP / detection frequency) as modeled in the simATAC
framework (Navidi et al., Genome Biology 2021): NZP = c0 + c1 \* mean +
c2 \* mean^2.

## Usage

``` r
calc_accessibility_sparsity_curve(data, poly_degree = 2)
```

## Arguments

- data:

  Count matrix (features/peaks x cells) or data.frame with 'mean' and
  'nzp'.

- poly_degree:

  Degree of polynomial (default: 2 for quadratic curve).

## Value

A list with estimated coefficients (c0, c1, c2), R-squared,
Spearman/Pearson correlations, and model fit summary.

## Examples

``` r
data(example_scatac, package = "scSimEval")
calc_accessibility_sparsity_curve(example_scatac$ref)
#> $c0
#> [1] -0.006022981
#> 
#> $c1
#> [1] 1.076984
#> 
#> $c2
#> [1] -0.388117
#> 
#> $r_squared
#> [1] 0.9692997
#> 
#> $spearman_cor
#> [1] 0.9781859
#> 
#> $pearson_cor
#> [1] 0.9783617
#> 
#> $poly_degree
#> [1] 2
#> 
#> $model
#> 
#> Call:
#> stats::lm(formula = y ~ stats::poly(x, degree = poly_degree, 
#>     raw = TRUE), data = df)
#> 
#> Coefficients:
#>                                       (Intercept)  
#>                                         -0.006023  
#> stats::poly(x, degree = poly_degree, raw = TRUE)1  
#>                                          1.076984  
#> stats::poly(x, degree = poly_degree, raw = TRUE)2  
#>                                         -0.388117  
#> 
#> 
```
