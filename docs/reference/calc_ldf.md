# Local Density Factor (LDF)

Computes the Local Density Estimate (LDE) and Local Density Factor (LDF)
using a Gaussian kernel over reachability distances in the k-NN
neighborhood, adapted from Latecki et al. and CellMixS.

## Usage

``` r
calc_ldf(coords, k = 15, h = 1, c = 1)
```

## Arguments

- coords:

  Matrix of coordinates / embeddings (cells x dimensions).

- k:

  Number of nearest neighbors. Default is 15.

- h:

  Bandwidth parameter for Gaussian kernel. Default is 1.

- c:

  Scaling constant for comparison of LDE to neighboring observations.
  Default is 1.

## Value

A named list:

- lde:

  Local density estimate for each cell

- ldf:

  Local density factor for each cell

## Examples

``` r
coords <- matrix(stats::rnorm(100), 50, 2)
batch <- factor(rep(c("B1", "B2"), length.out = 50))
calc_ldf(coords)
#> $lde
#>  [1] 0.12661222 0.13830896 0.05796056 0.10289912 0.05263022 0.12449512
#>  [7] 0.11091917 0.13999359 0.08745553 0.12778725 0.04956140 0.06103576
#> [13] 0.14071537 0.13059045 0.04639527 0.07602115 0.13575095 0.12174439
#> [19] 0.09051649 0.07389256 0.06159685 0.13453622 0.04653223 0.02548760
#> [25] 0.09433550 0.06076777 0.07987570 0.05178686 0.11564088 0.06102229
#> [31] 0.06129319 0.11082350 0.11460315 0.10695080 0.12913701 0.13290313
#> [37] 0.06293559 0.13888158 0.04084997 0.05852344 0.06336943 0.05348191
#> [43] 0.05302038 0.06825045 0.02132903 0.13695098 0.05316768 0.07448133
#> [49] 0.07721157 0.10865473
#> 
#> $ldf
#>  [1] 0.5009318 0.4801173 0.6164703 0.5191519 0.6554622 0.4958257 0.4991920
#>  [8] 0.4768767 0.5214887 0.4941199 0.6286867 0.6192561 0.4764290 0.4958125
#> [15] 0.6298302 0.5444494 0.4783714 0.4981606 0.5320438 0.5797362 0.5743400
#> [22] 0.4886118 0.6829581 0.7850851 0.5228667 0.5979899 0.5413838 0.5998397
#> [29] 0.5005631 0.5927448 0.6091571 0.4994222 0.4931281 0.4993083 0.4987990
#> [36] 0.4838091 0.5959644 0.4758096 0.6774077 0.6450886 0.5968440 0.6413921
#> [43] 0.6473359 0.5753667 0.7917231 0.4782941 0.6404251 0.5741198 0.5535384
#> [50] 0.5089279
#> 
```
