# Inverse Simpson Index (ISI) for Batch Mixing

Evaluates the diversity of batch labels within each cell's k-nearest
neighborhood using the Inverse Simpson Index (1 / sum(p_b^2)).
Implements a pure-R, lightweight formulation of LISI with optional
distance weighting.

## Usage

``` r
calc_isi(coords, batch_info, k = 30, weighted = TRUE)
```

## Arguments

- coords:

  Matrix of cell coordinates (cells x dimensions).

- batch_info:

  Factor or character vector of batch assignments.

- k:

  Number of nearest neighbors. Default is 30.

- weighted:

  Logical, whether to weight neighbor contributions by inverse distance.
  Default is TRUE.

## Value

A named list:

- isi_scores:

  Numeric vector of cell-level ISI values

- mean_isi:

  Mean ISI score (ranges from 1 to number of batches; higher indicates
  better mixing)

- median_isi:

  Median ISI score

## Examples

``` r
coords <- matrix(stats::rnorm(100), 50, 2)
batch <- factor(rep(c("B1", "B2"), length.out = 50))
calc_isi(coords, batch)
#> $isi_scores
#>  [1] 1.927025 1.938315 1.991476 1.940916 1.735261 1.884688 1.841849 1.928222
#>  [9] 1.891884 1.916513 1.952115 1.944333 1.819865 1.924779 1.996025 1.911742
#> [17] 1.915590 1.891906 1.857968 1.966656 1.824496 1.978098 1.899823 1.997787
#> [25] 1.938739 1.799605 1.990124 1.988413 1.991868 1.977942 1.960569 1.977269
#> [33] 1.820889 1.994421 1.969023 1.955415 1.958177 1.949180 1.974654 1.911693
#> [41] 1.984612 1.989831 1.868620 1.972950 1.963773 1.890808 1.954281 1.956124
#> [49] 1.996351 1.934521
#> 
#> $mean_isi
#> [1] 1.930944
#> 
#> $median_isi
#> [1] 1.946756
#> 
```
