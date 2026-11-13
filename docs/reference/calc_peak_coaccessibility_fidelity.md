# Evaluate Single-Cell ATAC Peak Co-Accessibility Fidelity (SCRIP)

Evaluates the preservation of chromatin peak-to-peak co-accessibility
correlation matrices (cis-regulatory interactions) between reference and
simulated scATAC-seq datasets.

## Usage

``` r
calc_peak_coaccessibility_fidelity(
  ref_atac,
  sim_atac,
  top_n_peaks = 300,
  method = c("spearman", "pearson")
)
```

## Arguments

- ref_atac:

  Reference scATAC-seq matrix (peaks x cells).

- sim_atac:

  Simulated scATAC-seq matrix (peaks x cells).

- top_n_peaks:

  Number of highest variance peaks to evaluate (default 300).

- method:

  Correlation method: "spearman" (default) or "pearson".

## Value

A list containing RV coefficient, matrix correlation, Frobenius
distance, and MAE.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_peak_coaccessibility_fidelity(example_scrna$ref, example_scrna$sim)
#> $rv_coefficient
#> [1] 0.5693799
#> 
#> $coaccessibility_pearson
#> [1] -0.01329664
#> 
#> $coaccessibility_spearman
#> [1] -0.008190968
#> 
#> $frobenius_distance
#> [1] 0.005381708
#> 
#> $coaccessibility_mae
#> [1] 0.1274649
#> 
```
