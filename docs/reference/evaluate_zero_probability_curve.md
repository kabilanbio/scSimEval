# Evaluate Zero-Probability Dropout Curve Fidelity (Splatter & ZINB-WaVE)

Compares empirical logistic dropout trends between reference and
simulated scRNA-seq datasets.

## Usage

``` r
evaluate_zero_probability_curve(ref_counts, sim_counts)
```

## Arguments

- ref_counts:

  Reference count matrix.

- sim_counts:

  Simulated count matrix.

## Value

Named numeric vector of curve parameter differences (`slope_error`,
`midpoint_error`, `r2_discrepancy`).

## Examples

``` r
data(example_scrna, package = "scSimEval")
evaluate_zero_probability_curve(example_scrna$ref, example_scrna$sim)
#>    slope_error midpoint_error r2_discrepancy 
#>      0.2560909      0.2624144      0.1646825 
```
