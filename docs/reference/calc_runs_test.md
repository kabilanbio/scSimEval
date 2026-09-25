# Calculate Wald-Wolfowitz Runs Test on Pooled Distributions

Integrated from countsimQC (Soneson & Robinson). Tests whether values
from the two datasets intermingle randomly when sorted. A significant
left-sided result indicates that identical values cluster together into
fewer runs than expected by chance, signaling distinct distributions.

## Usage

``` r
calc_runs_test(ref, sim, alternative = c("left.sided", "two.sided"))
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

- alternative:

  Alternative hypothesis: "left.sided" (default in countsimQC) or
  "two.sided".

## Value

A named list with the runs test statistic and p-value.
