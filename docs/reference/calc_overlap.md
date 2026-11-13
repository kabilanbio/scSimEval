# Calculate Kernel Density Overlap Value (OV)

Integrated from simpipe (Sun et al., Genome Biology 2021;
duohongrui/simpipe). Quantifies the shared area under the empirical
probability density curves of the reference and simulated distributions:
\$\$\mathrm{OV}(P, Q) = \int\_{-\infty}^{+\infty}
\min(\hat{f}\_{\mathrm{ref}}(x), \hat{f}\_{\mathrm{sim}}(x))\\dx\$\$
Uses
[`overlapping::overlap()`](https://rdrr.io/pkg/overlapping/man/overlap.html)
when available, with an exact normalized kernel density trapezoidal
integration fallback.

## Usage

``` r
calc_overlap(ref, sim, nbins = 1024)
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

- nbins:

  Number of equally spaced grid points for density integration (default
  1024, as in simpipe).

## Value

Numeric scalar in `[0, 1]` representing distributional overlap (1 =
identical distributions, 0 = disjoint; higher is better).

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_overlap(ref, sim)
#> [1] 0.8881093
```
