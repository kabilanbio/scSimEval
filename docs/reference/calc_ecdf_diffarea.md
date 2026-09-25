# Calculate Area Between Empirical Cumulative Distribution Functions (eCDFs)

Integrated from countsimQC (Soneson & Robinson). Measures the normalized
area between the eCDFs of reference and simulation across the shared
support.

## Usage

``` r
calc_ecdf_diffarea(ref, sim)
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

## Value

Normalized area between the two eCDFs.
