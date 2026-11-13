# Calculate Lineage Gene Expression Pattern Fidelity

Evaluates the preservation of continuous gene expression trends along
developmental trajectories by calculating the mean correlation of
expression profiles across pseudotime for top dynamic lineage genes.

## Usage

``` r
calc_lineage_expression_fidelity(
  ref_data,
  sim_data,
  ref_pseudotime = NULL,
  sim_pseudotime = NULL,
  n_genes = 50,
  n_bins = 20
)
```

## Arguments

- ref_data:

  Reference expression/count matrix.

- sim_data:

  Simulated expression/count matrix.

- ref_pseudotime:

  Optional precomputed reference pseudotime.

- sim_pseudotime:

  Optional precomputed simulated pseudotime.

- n_genes:

  Number of dynamic genes to evaluate (default 50).

- n_bins:

  Number of pseudotime trajectory bins (default 20).

## Value

Mean correlation of gene expression trends along lineages in \[0, 1\].

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_lineage_expression_fidelity(example_scrna$ref, example_scrna$sim)
#> [1] 0.1141288
```
