# Calculate Delta Variance (Pseudoreplication Bias Metric)

Measures the difference in gene expression variance between true
biological replicates and cell-shuffled pseudoreplicates (Squair et al.,
Nature Communications 2021). Positive delta variance indicates that
biological replicates have higher inter-individual variability than
randomly pooled pseudoreplicates.

## Usage

``` r
calc_delta_variance(
  counts,
  replicates,
  conditions,
  cell_types = NULL,
  n_perm = 3
)
```

## Arguments

- counts:

  Count matrix (features x cells).

- replicates:

  Factor of biological replicate / batch assignments for each cell.

- conditions:

  Factor of condition or cell-type labels.

- cell_types:

  Optional cell type vector to evaluate per cell type or overall.

- n_perm:

  Number of permutation shuffles to compute mean pseudoreplicate
  variance. Default is 3.

## Value

A named list:

- delta_variance:

  Numeric vector of gene-level delta variance values

- mean_delta_variance:

  Mean delta variance across all genes

- median_delta_variance:

  Median delta variance across all genes

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_delta_variance(example_scrna$ref,
                    replicates = example_scrna$batch_info,
                    conditions = example_scrna$cell_types)
#> $delta_variance
#>      Gene_01      Gene_02      Gene_03      Gene_04      Gene_05      Gene_06 
#>  -8805449.16  -7797415.07  24081491.72   3468783.90   6796240.15   -505415.72 
#>      Gene_07      Gene_08      Gene_09      Gene_10      Gene_11      Gene_12 
#>   9957176.85   -573202.00  33802655.33 -35716036.40   9079689.30 -31740032.93 
#>      Gene_13      Gene_14      Gene_15      Gene_16      Gene_17      Gene_18 
#>   5794475.71  -1874128.62   2619460.97 -13741501.18    629707.79  -2719745.01 
#>      Gene_19      Gene_20      Gene_21      Gene_22      Gene_23      Gene_24 
#>   2169671.29 -13633019.58   9908870.76   1710354.95   1553387.95  15139170.86 
#>      Gene_25      Gene_26      Gene_27      Gene_28      Gene_29      Gene_30 
#>  17741301.10   4994699.49 -20459354.80   2973031.82  11541643.36   5041622.22 
#>      Gene_31      Gene_32      Gene_33      Gene_34      Gene_35      Gene_36 
#>   5439507.51   4903711.25   1273624.73  -4282445.81   1198497.09   3200113.73 
#>      Gene_37      Gene_38      Gene_39      Gene_40      Gene_41      Gene_42 
#>    -65680.85  -1644783.56    101090.29   -138426.52 -19655238.49  -6857433.86 
#>      Gene_43      Gene_44      Gene_45      Gene_46      Gene_47      Gene_48 
#>   8320601.82   3671546.93   8332198.63 -15149545.82   1658506.61   2254000.39 
#>      Gene_49      Gene_50      Gene_51      Gene_52      Gene_53      Gene_54 
#> -27358712.06  12013883.94   7991292.66 -37497369.26  -5664444.26  16824197.57 
#>      Gene_55      Gene_56      Gene_57      Gene_58      Gene_59      Gene_60 
#>   3891178.79 -13662689.35 -43171380.89  -5417150.23    847173.15  16505430.65 
#> 
#> $mean_delta_variance
#> [1] -845010.2
#> 
#> $median_delta_variance
#> [1] 1605947
#> 
```
