# Perform Likelihood Ratio Test for Comparing Nested Single-Cell Simulation Models

Direct port and generalization of
[`scDesign3::perform_lrt`](https://rdrr.io/pkg/scDesign3/man/perform_lrt.html)
(Song et al., Nat Biotechnol 2024). Performs the likelihood ratio test
to compare two nested simulation models (e.g., cell-type/covariate model
vs intercept-only null model, or spline trajectory vs linear model).

## Usage

``` r
calc_likelihood_ratio_test(
  alter_model,
  null_model,
  df_alter = NULL,
  df_null = NULL
)
```

## Arguments

- alter_model:

  Alternative model (more complex) or list of alternative models per
  gene, or numeric log-likelihoods.

- null_model:

  Null model (simpler, strictly nested) or list of null models per gene,
  or numeric log-likelihoods.

- df_alter:

  Degrees of freedom for alternative model (used if models are numeric
  log-likelihoods).

- df_null:

  Degrees of freedom for null model (used if models are numeric
  log-likelihoods).

## Value

A data.frame containing:

- LogLik_alter:

  Log-likelihood under alternative model.

- LogLik_null:

  Log-likelihood under null model.

- df_alter:

  Degrees of freedom of alternative model.

- df_null:

  Degrees of freedom of null model.

- LR_statistic:

  Likelihood ratio statistic (-2 \* (LL_null - LL_alter)).

- delta_df:

  Difference in degrees of freedom (df_alter - df_null).

- p_value:

  P-value from chi-squared test with delta_df degrees of freedom.
