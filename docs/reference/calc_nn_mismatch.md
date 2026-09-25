# Calculate Nearest-Neighbor Label Mismatch Proportion

Integrated from countsimQC (Soneson & Robinson). Evaluates whether the
dataset label composition in k-NN neighborhoods departs significantly
from the overall global dataset proportion using Chi-squared tests.
Works for both 1D numeric vectors and 2D matrices.

## Usage

``` r
calc_nn_mismatch(ref, sim, k = NULL, subsample_size = 300)
```

## Arguments

- ref:

  Vector or matrix for reference dataset.

- sim:

  Vector or matrix for simulated dataset.

- k:

  Number of nearest neighbors (default max(5, 0.05 \* N)).

- subsample_size:

  Number of subsampled points to test (default 300).

## Value

Fraction of points with significant neighbor composition mismatch (p \<=
0.05).
