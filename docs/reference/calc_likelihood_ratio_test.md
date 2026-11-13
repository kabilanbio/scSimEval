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

## Examples

``` r
calc_likelihood_ratio_test(stats::rnorm(50), stats::rnorm(50))
#>    LogLik_alter LogLik_null df_alter df_null LR_statistic delta_df p_value
#> 1   -0.20135040 -0.97770728       NA      NA   1.55271375       NA      NA
#> 2   -0.28135682  0.62642572       NA      NA  -1.81556508       NA      NA
#> 3    0.21457297 -0.30603213       NA      NA   1.04121021       NA      NA
#> 4   -0.24654339 -0.64066865       NA      NA   0.78825051       NA      NA
#> 5   -1.06055986  0.26099615       NA      NA  -2.64311202       NA      NA
#> 6   -1.37503339  0.42976556       NA      NA  -3.60959790       NA      NA
#> 7   -0.87034897 -0.75851188       NA      NA  -0.22367419       NA      NA
#> 8   -1.26596173 -1.87715577       NA      NA   1.22238809       NA      NA
#> 9    0.49904436  0.49255531       NA      NA   0.01297811       NA      NA
#> 10  -0.13999982  0.81506228       NA      NA  -1.91012421       NA      NA
#> 11  -0.32542001  0.24751206       NA      NA  -1.14586414       NA      NA
#> 12   0.07240789 -1.12081233       NA      NA   2.38644044       NA      NA
#> 13   0.39198315  0.41252012       NA      NA  -0.04107393       NA      NA
#> 14  -1.35576887  0.29047640       NA      NA  -3.29249054       NA      NA
#> 15   1.23429759 -0.10609951       NA      NA   2.68079419       NA      NA
#> 16  -0.19216042  0.76938888       NA      NA  -1.92309859       NA      NA
#> 17   0.79862927  0.87080932       NA      NA  -0.14436011       NA      NA
#> 18   0.95737185 -1.03743347       NA      NA   3.98961065       NA      NA
#> 19   0.15294566 -2.06553215       NA      NA   4.43695562       NA      NA
#> 20   0.17516593 -1.53839812       NA      NA   3.42712810       NA      NA
#> 21  -0.12329344 -0.49097064       NA      NA   0.73535440       NA      NA
#> 22  -0.04801392  1.86999247       NA      NA  -3.83601278       NA      NA
#> 23   1.75629862 -0.66042978       NA      NA   4.83345681       NA      NA
#> 24  -1.65449233 -1.25667141       NA      NA  -0.79564183       NA      NA
#> 25   1.34977115  0.58866592       NA      NA   1.52221046       NA      NA
#> 26  -0.88479087  0.85867101       NA      NA  -3.48692377       NA      NA
#> 27  -0.54047500  0.60153036       NA      NA  -2.28401071       NA      NA
#> 28  -1.35811809 -1.86347368       NA      NA   1.01071118       NA      NA
#> 29   1.22424838 -1.02334124       NA      NA   4.49517925       NA      NA
#> 30   0.86658508  0.29245417       NA      NA   1.14826181       NA      NA
#> 31  -0.09836715 -1.79867894       NA      NA   3.40062357       NA      NA
#> 32  -1.24712396 -0.20772957       NA      NA  -2.07878879       NA      NA
#> 33   0.73815876 -0.76361207       NA      NA   3.00354166       NA      NA
#> 34   0.16388576  0.96839153       NA      NA  -1.60901154       NA      NA
#> 35  -0.10370943 -0.46694228       NA      NA   0.72646570       NA      NA
#> 36  -1.05310136 -1.12368349       NA      NA   0.14116426       NA      NA
#> 37  -0.43074005 -0.80793744       NA      NA   0.75439478       NA      NA
#> 38  -0.01307122 -0.22766684       NA      NA   0.42919124       NA      NA
#> 39  -0.89361693  0.32480076       NA      NA  -2.43683539       NA      NA
#> 40   0.04953601  1.71140423       NA      NA  -3.32373643       NA      NA
#> 41   0.01453970 -0.16695623       NA      NA   0.36299187       NA      NA
#> 42  -1.06905941 -1.50982461       NA      NA   0.88153040       NA      NA
#> 43   0.86849690  3.18362106       NA      NA  -4.63024834       NA      NA
#> 44   0.43383222  0.30994572       NA      NA   0.24777299       NA      NA
#> 45  -0.64841041 -0.08723892       NA      NA  -1.12234299       NA      NA
#> 46   0.28569106  2.37648576       NA      NA  -4.18158940       NA      NA
#> 47  -1.29975719  0.41118264       NA      NA  -3.42187967       NA      NA
#> 48   0.51602401  0.22042373       NA      NA   0.59120057       NA      NA
#> 49   1.12424620 -0.19882059       NA      NA   2.64613356       NA      NA
#> 50  -1.26035271 -0.40073659       NA      NA  -1.71923224       NA      NA
```
