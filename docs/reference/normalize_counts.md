# Fast Count Normalization (Log2-CPM or TF-IDF)

Converts raw single-cell count matrices into log-normalized CPM
(`"lognorm"`) for scRNA-seq or TF-IDF (`"tfidf"`) for scATAC-seq,
preserving `dgCMatrix` sparsity.

## Usage

``` r
normalize_counts(
  mat,
  scale_factor = 1e+06,
  log_base = 2,
  pseudo_count = 1,
  method = c("lognorm", "tfidf")
)
```

## Arguments

- mat:

  A count matrix (genes/features x cells), `dgCMatrix`,
  `SingleCellExperiment`, or `Seurat` object.

- scale_factor:

  Numeric scaling factor for counts per cell. Default is 1e6 (CPM) for
  `"lognorm"`, or 1e4 for `"tfidf"`.

- log_base:

  Numeric base for the logarithm in `"lognorm"` mode: 2 (default, for
  log2(CPM+1)) or natural log `exp(1)`.

- pseudo_count:

  Numeric pseudo-count added before logarithm in `"lognorm"` mode.
  Default is 1.

- method:

  Normalization method: `"lognorm"` (default) or `"tfidf"`
  (Signac/simPIC scATAC-seq standard).

## Value

Normalized matrix of transformed values.

## Examples

``` r
data(example_scrna, package = "scSimEval")
norm_mat <- normalize_counts(example_scrna$ref)
```
