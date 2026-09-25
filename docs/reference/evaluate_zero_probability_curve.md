# Evaluate Zero-Probability Dropout Curve Fidelity (Splatter & ZINB-WaVE)

Compares empirical dropout trends between reference and simulated
datasets.

## Usage

``` r
evaluate_zero_probability_curve(ref_counts, sim_counts)
```

## Arguments

- ref_counts:

  Reference count matrix.

- sim_counts:

  Simulated count matrix.

## Value

Named numeric vector of curve parameter differences.
