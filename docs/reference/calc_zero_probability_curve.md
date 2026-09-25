# Fit Zero-Probability Dropout Curve (Splatter & ZINB-WaVE)

Fits an empirical logistic dropout curve: logit(P(Y=0)) = beta_0 +
beta_1 \* log(mu) to model the dropout relationship with mean
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
