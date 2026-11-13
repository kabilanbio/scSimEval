# Fast BLAS-Accelerated Cross-Distance Squared Matrix

Computes \|\|X_i - Y_j\|\|^2 using BLAS matrix product: \|\|x\|\|^2 +
\|\|y\|\|^2 - 2 \<x, y\>.

## Usage

``` r
fast_cross_distance_sq(X, Y)
```

## Arguments

- X:

  Matrix (m x p).

- Y:

  Matrix (n x p).

## Value

Matrix of squared Euclidean distances (m x n).
