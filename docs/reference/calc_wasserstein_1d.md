# Calculate 1D Wasserstein Metric / Earth Mover's Distance (WS) Integrated from HelenaLC/simulation-comparison.

Calculate 1D Wasserstein Metric / Earth Mover's Distance (WS) Integrated
from HelenaLC/simulation-comparison.

## Usage

``` r
calc_wasserstein_1d(ref, sim, p = 1)
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

- p:

  Power of the Wasserstein metric (default 1 for standard earth mover's
  distance).

## Value

The 1D Wasserstein distance between the two empirical distributions.
