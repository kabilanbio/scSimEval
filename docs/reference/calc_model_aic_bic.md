# Compute Model Information Criteria (AIC and BIC)

Evaluates the statistical trade-off between model fit and parameter
complexity for single-cell generative models as benchmarked in scDesign3
(Song et al., Nat Biotechnol 2024).

## Usage

``` r
calc_model_aic_bic(loglik, n_params, n_obs)
```

## Arguments

- loglik:

  Numeric vector or scalar of log-likelihood values.

- n_params:

  Numeric vector or scalar of number of estimated parameters (degrees of
  freedom).

- n_obs:

  Total number of independent observations (cells or cell-gene pairs).

## Value

A named vector or data.frame containing loglik, n_params, n_obs, AIC,
and BIC.
