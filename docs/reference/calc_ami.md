# Calculate Adjusted Mutual Information (AMI)

Computes the Adjusted Mutual Information (AMI) between two clusterings
or cell-type label assignments, correcting for chance expected mutual
information under the hypergeometric model with fixed marginals (Vinh et
al., JMLR 2010), as used in scMoMtF (Lan et al., PLOS Comput Biol 2024).

## Usage

``` r
calc_ami(
  pred,
  truth,
  average_method = c("arithmetic", "max", "min", "geometric")
)
```

## Arguments

- pred:

  Vector of predicted cluster labels.

- truth:

  Vector of ground truth cell type labels.

- average_method:

  Method to normalize mutual information: "arithmetic" (default,
  scikit-learn standard), "max", "min", or "geometric".

## Value

AMI score in \[0, 1\] (adjusted for chance).

## References

Vinh, N. X., Epps, J., & Bailey, J. (2010). Information theoretic
measures for clusterings comparison: Variants, properties, normalization
and correction for chance. Journal of Machine Learning Research, 11,
2837-2854.

Lan, W., Ling, T., Chen, Q. et al. scMoMtF: An interpretable multitask
learning framework for single-cell multi-omics data analysis. PLOS
Comput Biol 20(12): e1012679 (2024).
