# Detect Differential Variability (DV) Genes (Bartlett Test)

Detect Differential Variability (DV) Genes (Bartlett Test)

## Usage

``` r
calc_signal_dv(exprs_mat, cell_types)
```

## Arguments

- exprs_mat:

  Log-normalized matrix (genes x cells).

- cell_types:

  Factor of 2 cell types.

## Value

Vector of adjusted p-values.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_signal_dv(example_scrna$ref, example_scrna$cell_types)
#>  [1] 7.656429e-02 2.756860e-02 1.754115e-05 2.258443e-01 2.010435e-01
#>  [6] 8.016641e-01 4.553276e-01 9.328213e-01 7.172070e-01 7.656429e-02
#> [11] 2.573023e-02 2.708077e-02 1.337701e-07 3.910768e-02 7.025884e-03
#> [16] 9.534988e-01 9.534988e-01 9.534988e-01 9.534988e-01 2.935425e-02
#> [21] 7.522318e-02 2.292247e-01 7.656429e-02 3.644857e-01 9.359630e-01
#> [26] 8.240211e-01 7.172070e-01 2.665707e-01 8.240211e-01 8.240211e-01
#> [31] 3.644857e-01 7.656429e-02 8.639075e-01 3.850341e-02 8.803781e-03
#> [36] 4.759268e-01 4.344661e-01 8.803781e-03 2.979180e-02 8.240211e-01
#> [41] 9.909546e-01 1.636175e-01 4.759268e-01 9.433383e-01 1.636175e-01
#> [46] 3.850341e-02 8.600443e-01 8.240211e-01 8.240211e-01 4.238203e-01
#> [51] 6.828151e-01 6.028387e-01 4.344661e-01 6.828151e-01 1.322160e-03
#> [56] 7.858415e-01 3.905416e-01 4.759268e-01 1.306924e-02 2.258443e-01
```
