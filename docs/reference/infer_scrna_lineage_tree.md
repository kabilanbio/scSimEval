# Automatically Infer Lineage Tree (hclust) from scRNA-seq Counts & Cell Types

Builds a hierarchical clustering lineage tree between cell-type
centroids from scRNA-seq expression data.

## Usage

``` r
infer_scrna_lineage_tree(data, cell_types, method = "ward.D2")
```

## Arguments

- data:

  Expression or count matrix (genes x cells).

- cell_types:

  Factor or vector of cell-type annotations for each cell.

- method:

  Linkage method for hierarchical clustering (default "ward.D2").

## Value

An `hclust` tree object representing the cell-type lineage tree.
