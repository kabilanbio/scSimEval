# Marginal Model Goodness of Fit and Information Criteria for Single-Cell Simulators

Evaluates gene-wise and aggregate marginal goodness-of-fit (Poisson,
Negative Binomial, or Gaussian) for simulated or fitted single-cell
count matrices, computing total and mean AIC and BIC across genes
(scDesign3; Song et al., 2024).

## Usage

``` r
calc_marginal_aic_bic(
  counts,
  fitted_means,
  dispersions = NULL,
  distribution = c("poisson", "nb", "gaussian"),
  n_params_per_gene = NULL
)
```

## Arguments

- counts:

  Matrix or data.frame of observed or simulated counts (genes x cells).

- fitted_means:

  Matrix or data.frame of model fitted means or expectations (genes x
  cells).

- dispersions:

  Optional vector of gene-level dispersion parameters for Negative
  Binomial model. If NULL, estimated via method-of-moments.

- distribution:

  Parametric distribution: "poisson", "nb" (Negative Binomial), or
  "gaussian".

- n_params_per_gene:

  Number of estimated parameters per gene. If NULL, defaults to 1 for
  Poisson, 2 for NB, and 2 for Gaussian.

## Value

A list containing:

- total_loglik:

  Sum of marginal log-likelihoods across all genes and cells.

- total_aic:

  Aggregate AIC across all genes.

- total_bic:

  Aggregate BIC across all genes.

- mean_gene_aic:

  Mean AIC per gene.

- mean_gene_bic:

  Mean BIC per gene.

- gene_loglik:

  Vector of log-likelihoods per gene.

- gene_aic:

  Vector of AIC per gene.

- gene_bic:

  Vector of BIC per gene.
