# Evaluate the 5 SimBench Biological Signal Proportions

Compares proportions of genes exhibiting DE, DV, DD, DP, and BD between
empirical reference and simulated datasets.

## Usage

``` r
evaluate_simbench_signals(
  ref_mat,
  sim_mat,
  ref_celltypes,
  sim_celltypes,
  p_sig = 0.05,
  bi_cutoff = 0.3
)
```

## Arguments

- ref_mat:

  Reference matrix.

- sim_mat:

  Simulated matrix.

- ref_celltypes:

  Reference cell type labels (subsets to top 2 abundant types).

- sim_celltypes:

  Simulated cell type labels.

- p_sig:

  Significance cutoff (default 0.05).

- bi_cutoff:

  Bimodal index cutoff (default 0.3).

## Value

A tidy data.frame comparing biological signal proportions.

## Examples

``` r
data(example_scrna, package = "scSimEval")
evaluate_simbench_signals(example_scrna$ref, example_scrna$sim,
                          ref_celltypes = example_scrna$cell_types,
                          sim_celltypes = example_scrna$cell_types)
#>               Signal_Type Reference_Prop Simulation_Prop Absolute_Error
#> 1         DE (Mean Shift)     0.00000000      0.00000000     0.00000000
#> 2        DV (Variability)     0.06666667      0.01666667     0.05000000
#> 3       DD (Distribution)     0.00000000      0.03333333     0.03333333
#> 4 DP (Proportion/Dropout)     0.00000000      0.00000000     0.00000000
#> 5      BD (Bimodal Index)     0.25000000      0.20000000     0.05000000
```
