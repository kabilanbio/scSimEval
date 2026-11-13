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

## Examples

``` r
calc_marginal_aic_bic(stats::rnorm(50), stats::rnorm(50))
#> $total_loglik
#> [1] -85.59207
#> 
#> $total_aic
#> [1] 271.1841
#> 
#> $total_bic
#> [1] 171.1841
#> 
#> $mean_gene_aic
#> [1] 5.423683
#> 
#> $mean_gene_bic
#> [1] 3.423683
#> 
#> $gene_loglik
#>  [1]  -0.7254615 -12.5478097  -8.5869897  -0.9867235  35.5428058  -8.4363095
#>  [7]   3.6279682 -13.3151493  -1.3637766  -2.5792225 -23.6198578 -23.6208842
#> [13]  -8.2002859  16.1762422  -0.6105219  -3.2961656  -0.1402087 -16.2981509
#> [19]  -0.2224286 -32.4619893   3.4769188  12.9897768  -1.2165843 -12.4133247
#> [25]  -1.1898062   4.4802504  -0.6407442  -1.4212662  -2.7923962  -0.3951945
#> [31]  53.8950137  -0.3659452  -4.0011372  -1.6225933  36.6961798  -0.6972320
#> [37] -21.7789427  -8.3486581  -1.9355572  32.8709180  -0.9837802  -0.7037743
#> [43]  -3.8463705   7.3430071  -3.1693101  10.5943243  -1.9616582  -3.2300200
#> [49] -22.3687291 -51.1905179
#> 
#> $gene_aic
#>  [1]    3.450923   27.095619   19.173979    3.973447  -69.085612   18.872619
#>  [7]   -5.255936   28.630299    4.727553    7.158445   49.239716   49.241768
#> [13]   18.400572  -30.352484    3.221044    8.592331    2.280417   34.596302
#> [19]    2.444857   66.923979   -4.953838  -23.979554    4.433169   26.826649
#> [25]    4.379612   -6.960501    3.281488    4.842532    7.584792    2.790389
#> [31] -105.790027    2.731890   10.002274    5.245187  -71.392360    3.394464
#> [37]   45.557885   18.697316    5.871114  -63.741836    3.967560    3.407549
#> [43]    9.692741  -12.686014    8.338620  -19.188649    5.923316    8.460040
#> [49]   46.737458  104.381036
#> 
#> $gene_bic
#>  [1]    1.4509230   25.0956193   17.1739794    1.9734469  -71.0856116
#>  [6]   16.8726190   -7.2559365   26.6302986    2.7275532    5.1584450
#> [11]   47.2397156   47.2417685   16.4005718  -32.3524844    1.2210439
#> [16]    6.5923313    0.2804173   32.5963018    0.4448572   64.9239786
#> [21]   -6.9538376  -25.9795536    2.4331686   24.8266494    2.3796125
#> [26]   -8.9605008    1.2814884    2.8425324    5.5847923    0.7903891
#> [31] -107.7900274    0.7318905    8.0022743    3.2451867  -73.3923597
#> [36]    1.3944639   43.5578853   16.6973162    3.8711144  -65.7418360
#> [41]    1.9675605    1.4075487    7.6927410  -14.6860142    6.3386203
#> [46]  -21.1886487    3.9233165    6.4600400   44.7374581  102.3810357
#> 
```
