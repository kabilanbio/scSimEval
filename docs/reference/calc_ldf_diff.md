# Local Density Differences (ldfDiff)

Quantifies cell-specific changes in the Local Density Factor (LDF)
before and after data integration or between reference and simulated
datasets.

## Usage

``` r
calc_ldf_diff(coords_pre, coords_post, k = 15, h = 1, c = 1)
```

## Arguments

- coords_pre:

  Matrix of coordinates before integration / reference (cells x
  dimensions).

- coords_post:

  Matrix of coordinates after integration / simulated (cells x
  dimensions).

- k:

  Number of nearest neighbors. Default is 15.

- h:

  Bandwidth parameter for Gaussian kernel. Default is 1.

- c:

  Scaling constant. Default is 1.

## Value

A named list:

- diff:

  Numeric vector of absolute LDF differences per cell

- mean_ldf_diff:

  Mean absolute LDF difference (lower indicates better structure
  preservation)

- median_ldf_diff:

  Median absolute LDF difference

- ldf_pre:

  LDF in pre-integration / reference space

- ldf_post:

  LDF in post-integration / simulated space

## Examples

``` r
coords <- matrix(stats::rnorm(100), 50, 2)
batch <- factor(rep(c("B1", "B2"), length.out = 50))
calc_ldf_diff(coords, coords)
#> $diff
#>  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> [39] 0 0 0 0 0 0 0 0 0 0 0 0
#> 
#> $mean_ldf_diff
#> [1] 0
#> 
#> $median_ldf_diff
#> [1] 0
#> 
#> $ldf_pre
#>  [1] 0.5584664 0.4882385 0.5484005 0.4851901 0.5130473 0.7360220 0.6200072
#>  [8] 0.4905670 0.5020313 0.5244704 0.6245361 0.5425613 0.5719693 0.5061452
#> [15] 0.6459696 0.5283880 0.7756749 0.5274256 0.5254505 0.4977375 0.4952705
#> [22] 0.6689918 0.5628271 0.4854934 0.5225224 0.5798444 0.4959112 0.5493756
#> [29] 0.5380227 0.4843981 0.5780511 0.5991865 0.5305217 0.5335420 0.4880888
#> [36] 0.4808848 0.4941464 0.5687223 0.4908002 0.5774380 0.5488743 0.5027382
#> [43] 0.6712591 0.4853762 0.5124770 0.4851068 0.5302969 0.4914481 0.5390119
#> [50] 0.5912593
#> 
#> $ldf_post
#>  [1] 0.5584664 0.4882385 0.5484005 0.4851901 0.5130473 0.7360220 0.6200072
#>  [8] 0.4905670 0.5020313 0.5244704 0.6245361 0.5425613 0.5719693 0.5061452
#> [15] 0.6459696 0.5283880 0.7756749 0.5274256 0.5254505 0.4977375 0.4952705
#> [22] 0.6689918 0.5628271 0.4854934 0.5225224 0.5798444 0.4959112 0.5493756
#> [29] 0.5380227 0.4843981 0.5780511 0.5991865 0.5305217 0.5335420 0.4880888
#> [36] 0.4808848 0.4941464 0.5687223 0.4908002 0.5774380 0.5488743 0.5027382
#> [43] 0.6712591 0.4853762 0.5124770 0.4851068 0.5302969 0.4914481 0.5390119
#> [50] 0.5912593
#> 
```
