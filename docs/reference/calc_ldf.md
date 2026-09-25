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
