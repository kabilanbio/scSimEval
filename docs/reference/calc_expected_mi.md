# Expected Mutual Information for Cluster Comparison

Computes the exact expectation of mutual information under the
generalized hypergeometric model of randomness with fixed marginals
(Vinh et al., 2010).

## Usage

``` r
calc_expected_mi(a, b, N)
```

## Arguments

- a:

  Integer vector of row sums (marginal cluster sizes of partition 1).

- b:

  Integer vector of column sums (marginal cluster sizes of partition 2).

- N:

  Total number of items / cells.

## Value

Expected mutual information in nats.
