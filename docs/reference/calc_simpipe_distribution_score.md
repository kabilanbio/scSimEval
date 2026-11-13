# Calculate Simpipe P-value Uniformity Distribution Score (Approach 2)

Evaluates the reliability of candidate DEGs returned by simulation
methods according to the Simpipe benchmarking framework (Duo et al.,
Genome Biology 2024).

## Usage

``` r
calc_simpipe_distribution_score(
  sim_data,
  sim_celltypes,
  candidate_degs = NULL,
  n_bins = 10,
  alpha = 0.05,
  min_cells = 5
)
```

## Arguments

- sim_data:

  Simulated count matrix (genes x cells).

- sim_celltypes:

  Factor or character vector of simulated cell groups.

- candidate_degs:

  Optional list of character vectors of candidate DEGs per pair to
  remove. If NULL, edgeRQLFDetRate is first run on sim_data to identify
  candidate DEGs for removal.

- n_bins:

  Number of histogram bins for the Pearson chi-squared goodness-of-fit
  test (default 10).

- alpha:

  Significance threshold for accepting the uniform null hypothesis
  (default 0.05). If test p-value \> alpha, null is accepted (S = 1).

- min_cells:

  Minimum number of cells per group (default 5).

## Value

A list containing:

- distribution_score:

  Overall Simpipe distribution score (0 to 1)

- pairwise_results:

  Data frame with per-pair details: Group1, Group2, N_removed,
  N_remaining, Chisq_Stat, Chisq_Pval, Accepted_Null, S_i

- pvalues_remaining:

  Named list of P-value vectors for the remaining genes per pair

## Details

The hypothesis is that if there are no DEGs between two sample groups,
the false positive rate of detected DEGs should be approximately 0.05
and P-values should follow a uniform distribution \\U(0, 1)\\. Candidate
DEGs returned by the simulation method are directly removed from the
simulated gene expression matrix, and differential expression analysis
is performed on the remaining genes using `edgeRQLFDetRate`.

The Pearson chi-squared goodness-of-fit test is performed on the
remaining P-values under the null hypothesis of uniformity across 10
equal bins. If the null hypothesis is accepted (\\P \> 0.05\\), \\S_i =
1\\; otherwise \\S_i = 0\\.

For 2 groups (\\n = 2\\): \$\$\text{distribution score} = 0 \text{ or }
1\$\$ For 3 or more groups (\\n \> 2\\): \$\$\text{distribution score} =
\frac{1}{N} \sum\_{i=1}^N S_i\$\$ where \\N = \binom{n}{2}\\ is the
number of paired groups.

## Examples

``` r
data(example_scrna, package = "scSimEval")
res <- calc_simpipe_distribution_score(
  sim_data = example_scrna$sim,
  sim_celltypes = example_scrna$cell_types
)
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
res$distribution_score
#> [1] 1
```
