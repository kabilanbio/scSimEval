# Calculate Between-Dataset Silhouette Width

Integrated from countsimQC (Soneson & Robinson). Treats dataset identity
(reference vs. simulation) as cluster labels to test if the two datasets
separate into distinct clusters or remain well-mixed.

## Usage

``` r
calc_between_dataset_silhouette(ref, sim, subsample_size = 300)
```

## Arguments

- ref:

  Vector or matrix for reference dataset.

- sim:

  Vector or matrix for simulated dataset.

- subsample_size:

  Number of subsampled points to evaluate (default 300).

## Value

A named list with global and local between-dataset silhouette widths.
