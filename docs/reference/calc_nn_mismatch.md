# Calculate Nearest-Neighbor Composition / Local Neighborhood Similarity (NN Mismatch)

Integrated from countsimQC (Soneson & Robinson). Evaluates local
neighborhood composition by pooling reference and simulated observations
(1D vector or 2D matrix), finding the \\k\\ nearest neighbors for each
observation, and performing a Pearson Chi-squared goodness-of-fit test
comparing the local neighborhood dataset proportions against the global
dataset proportions. Returns the fraction of observations with \\p \le
0.05\\.

## Usage

``` r
calc_nn_mismatch(ref, sim, k = NULL, subsample_size = 300)
```

## Arguments

- ref:

  Numeric vector (1D) or 2-column matrix/data.frame (2D) of reference
  values.

- sim:

  Numeric vector (1D) or 2-column matrix/data.frame (2D) of simulated
  values.

- k:

  Integer number of nearest neighbors. If NULL, defaults to
  `max(5, round(0.05 * N))`.

- subsample_size:

  Maximum number of observations per dataset to subsample for speed
  (default 300).

## Value

Numeric scalar in `[0, 1]` representing the fraction of local
neighborhoods with significant dataset composition mismatch (0 =
well-mixed, 1 = completely separated).

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_nn_mismatch(ref, sim)
#> [1] 0.07
```
