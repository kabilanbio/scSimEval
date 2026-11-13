# Perform TF-IDF Normalization for scATAC-seq Peak Matrices

Computes Term Frequency–Inverse Document Frequency (TF-IDF)
normalization following the Signac (`RunTFIDF(..., method = 1)`) and
simPIC (Chugh et al., 2024) standard for scATAC-seq chromatin
accessibility matrices: \$\$\text{TF-IDF}\_{i,j} = \log\left(1 +
\frac{X\_{i,j}}{N_j} \cdot \frac{N\_{\text{cells}}}{\sum_k X\_{i,k}}
\cdot \text{scale\\factor}\right)\$\$

## Usage

``` r
tfidf_transform(mat, scale_factor = 10000)
```

## Arguments

- mat:

  A peak-by-cell count matrix, sparse `dgCMatrix`,
  `SingleCellExperiment`, or `Seurat` object.

- scale_factor:

  Numeric scaling factor applied prior to `log1p`. Default is `1e4`
  (10,000).

## Value

TF-IDF transformed matrix (preserving `dgCMatrix` sparsity when sparse
input is supplied).

## Examples

``` r
set.seed(42)
m <- matrix(rpois(100, 0.5), nrow = 10)
tfidf_mat <- tfidf_transform(m)
```
