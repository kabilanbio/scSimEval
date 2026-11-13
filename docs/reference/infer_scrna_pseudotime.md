# Automatically Infer Pseudotime Trajectory from scRNA-seq Counts

Estimates cell differentiation pseudotime directly from single-cell
expression counts using diffusion/principal curve projection along the
first principal component of top highly variable genes.

## Usage

``` r
infer_scrna_pseudotime(data, n_top = 500)
```

## Arguments

- data:

  Expression or count matrix (genes x cells), or SingleCellExperiment
  object.

- n_top:

  Number of top variable genes to use (default 500).

## Value

Numeric vector of pseudotime values in \[0, 1\] for each cell.

## Examples

``` r
data(example_scrna, package = "scSimEval")
infer_scrna_pseudotime(example_scrna$ref)
#>    Cell_01    Cell_02    Cell_03    Cell_04    Cell_05    Cell_06    Cell_07 
#> 0.41759739 0.45901927 0.49637168 0.77353790 0.34387925 0.84781040 0.52877867 
#>    Cell_08    Cell_09    Cell_10    Cell_11    Cell_12    Cell_13    Cell_14 
#> 0.54790473 0.42466216 0.39397755 0.34431179 0.69547915 0.33091467 0.38385344 
#>    Cell_15    Cell_16    Cell_17    Cell_18    Cell_19    Cell_20    Cell_21 
#> 0.67072136 0.23467611 0.60165132 0.44670669 0.85094013 0.38061329 0.28963526 
#>    Cell_22    Cell_23    Cell_24    Cell_25    Cell_26    Cell_27    Cell_28 
#> 0.36523250 0.54499857 0.81496222 0.35765557 0.49405864 0.29677919 0.66358160 
#>    Cell_29    Cell_30    Cell_31    Cell_32    Cell_33    Cell_34    Cell_35 
#> 0.35807340 0.42839411 0.29153147 0.84343713 0.57763568 0.80586828 0.61264403 
#>    Cell_36    Cell_37    Cell_38    Cell_39    Cell_40    Cell_41    Cell_42 
#> 0.75374225 0.40209081 0.72121703 0.58482378 0.70780559 0.87851194 1.00000000 
#>    Cell_43    Cell_44    Cell_45    Cell_46    Cell_47    Cell_48    Cell_49 
#> 0.01599112 0.66400402 0.82871228 0.33085210 0.57096273 0.59476202 0.60315576 
#>    Cell_50    Cell_51    Cell_52    Cell_53    Cell_54    Cell_55    Cell_56 
#> 0.41207545 0.72506061 0.54074057 0.05658093 0.40222307 0.00000000 0.48314164 
#>    Cell_57    Cell_58    Cell_59    Cell_60    Cell_61    Cell_62    Cell_63 
#> 0.79510404 0.54572867 0.78919845 0.78704550 0.67376888 0.70963373 0.47499017 
#>    Cell_64    Cell_65    Cell_66    Cell_67    Cell_68    Cell_69    Cell_70 
#> 0.31516190 0.16160572 0.28523016 0.54812465 0.18068369 0.53231035 0.58105737 
#>    Cell_71    Cell_72    Cell_73    Cell_74    Cell_75    Cell_76    Cell_77 
#> 0.82357821 0.59226123 0.47771115 0.48496802 0.23412951 0.69182896 0.58624743 
#>    Cell_78    Cell_79    Cell_80 
#> 0.38078131 0.93076680 0.48873806 
```
