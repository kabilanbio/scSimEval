# Calculate Overall Clustering-Like Separation of Datasets via Silhouette-Based Measures

Integrated from countsimQC (Soneson & Robinson). Evaluates overall
separation between reference and simulated distributions by treating
dataset identity (reference vs. simulated) as cluster labels on the
pooled observations and computing both the global average silhouette
width and the local neighborhood silhouette width (based on the closest
\\k = \max(5, \text{round}(0.05 N))\\ points in each dataset).

## Usage

``` r
calc_between_dataset_silhouette(ref, sim, subsample_size = 300)
```

## Arguments

- ref:

  Numeric vector (1D) or 2-column matrix/data.frame (2D) of reference
  values.

- sim:

  Numeric vector (1D) or 2-column matrix/data.frame (2D) of simulated
  values.

- subsample_size:

  Maximum number of observations per dataset to subsample for speed
  (default 300).

## Value

A named list with `between_dataset_silh_global` and
`between_dataset_silh_local`.

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_between_dataset_silhouette(ref, sim)
#> $between_dataset_silh_global
#> [1] -0.01126145
#> 
#> $between_dataset_silh_local
#> [1] -0.1681511
#> 
```
