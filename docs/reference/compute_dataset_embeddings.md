# Compute Low-Dimensional Embeddings for Reference and Simulated Single-Cell Datasets

Projects empirical biological reference and simulated count matrices
into low-dimensional representations (UMAP, t-SNE, or PCA) for
side-by-side and multi-method topological benchmarking. Automatically
performs library-size normalization, variance-based feature selection,
principal component analysis, and optional unsupervised clustering.

## Usage

``` r
compute_dataset_embeddings(
  reference,
  simulated,
  reduction = c("umap", "tsne", "pca"),
  n_pcs = 30,
  perplexity = 30,
  n_neighbors = 15,
  min_dist = 0.3,
  seed = 42,
  cell_types = NULL,
  batch = NULL
)
```

## Arguments

- reference:

  Biological reference count matrix (features x cells), dgCMatrix,
  SingleCellExperiment, or Seurat object.

- simulated:

  A single simulated count matrix or a named list of simulated count
  matrices (e.g. `list("Splatter" = sim1, "scDesign3" = sim2)`).

- reduction:

  Character string specifying the dimensionality reduction method:
  `"umap"` (default), `"tsne"`, or `"pca"`.

- n_pcs:

  Integer specifying the number of principal components to calculate
  (default: 30).

- perplexity:

  Numeric perplexity for t-SNE (default: 30; automatically adapted for
  small sample sizes).

- n_neighbors:

  Integer number of nearest neighbors for UMAP (default: 15).

- min_dist:

  Numeric minimum distance parameter for UMAP (default: 0.3).

- seed:

  Random seed for reproducibility (default: 42).

- cell_types:

  Optional factor or character vector of cell type annotations for
  reference cells (or named list if per-dataset).

- batch:

  Optional factor or character vector of batch annotations.

## Value

A tidy `data.frame` containing cell coordinates (`Dim1`, `Dim2`),
`Dataset` name, `Dataset_Type` ("Reference" vs "Simulated"),
`Cell_Type`, `Cluster`, `Library_Size`, `Detected_Features`, and
`Batch`.

## Examples

``` r
data(example_scrna, package = "scSimEval")
emb <- compute_dataset_embeddings(
  reference = example_scrna$ref,
  simulated = list("Splatter" = example_scrna$sim),
  reduction = "umap",
  cell_types = example_scrna$cell_types
)
#> Warning: You're computing too large a percentage of total singular values, use a standard svd instead.
#> Warning: You're computing too large a percentage of total singular values, use a standard svd instead.
head(emb)
#>   Cell_ID       Dim1       Dim2   Dataset      Role Dataset_Type Cell_Type
#> 1 Cell_01  2.1842285  0.3066108 Reference Reference    Reference     TypeA
#> 2 Cell_02 -0.6592609 -1.5866539 Reference Reference    Reference     TypeA
#> 3 Cell_03  1.4652766  1.5027129 Reference Reference    Reference     TypeA
#> 4 Cell_04 -1.8158411 -0.5365884 Reference Reference    Reference     TypeA
#> 5 Cell_05 -0.3859525  0.9999394 Reference Reference    Reference     TypeA
#> 6 Cell_06 -0.4865342 -0.7764473 Reference Reference    Reference     TypeA
#>     Cluster Library_Size Detected_Features   Batch
#> 1 Cluster_2          225                47 Batch 1
#> 2 Cluster_2          208                49 Batch 1
#> 3 Cluster_1          262                48 Batch 1
#> 4 Cluster_2          227                47 Batch 1
#> 5 Cluster_1          255                51 Batch 1
#> 6 Cluster_1          256                46 Batch 1
```
