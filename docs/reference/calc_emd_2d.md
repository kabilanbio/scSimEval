# Calculate 2D Earth Mover's Distance (2D EMD) Integrated from HelenaLC/simulation-comparison.

Computes 2D bivariate density via MASS::kde2d over a shared bounding
box, then measures optimal transport distance via emdist::emd2d.

## Usage

``` r
calc_emd_2d(ref_mat, sim_mat, n = 25)
```

## Arguments

- ref_mat:

  2-column numeric matrix for reference.

- sim_mat:

  2-column numeric matrix for simulation.

- n:

  Grid resolution for 2D density estimation. Default is 25.

## Value

Normalized 2D Earth Mover's Distance.
