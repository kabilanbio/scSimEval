# Calculate 1D Kernel Density Estimation Integrated Squared Error (KDE-ISE)

Integrated from SimBench (SydneyBioX/SimBench). Measures the Integrated
Squared Error (ISE) between the empirical kernel density estimates of
reference and simulated properties across the shared support: \$\$ISE =
\int (\hat{f}\_{ref}(x) - \hat{f}\_{sim}(x))^2 dx\$\$

## Usage

``` r
calc_kde_ise(ref, sim, n_points = 512)
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

- n_points:

  Number of points for density grid evaluation (default 512).

## Value

Integrated Squared Error (non-negative numeric value, 0 indicates
identical density shapes).

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_kde_ise(ref, sim)
#> [1] 0.01216429
```
