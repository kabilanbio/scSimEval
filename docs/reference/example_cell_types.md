# Example Cell-Type Annotations Vector

A named factor vector containing discrete cell-type labels ("TypeA",
"TypeB") for 80 single cells, matching the column names of
[`example_scrna`](https://kabilanbio.github.io/scSimEval/reference/example_scrna.md)
and
[`example_scatac`](https://kabilanbio.github.io/scSimEval/reference/example_scatac.md).
Useful as a reference template for formatting cell-type metadata inputs
in `scSimEval`.

## Usage

``` r
data(example_cell_types)
```

## Format

A named factor with 80 elements across 2 levels: "TypeA" (40 cells) and
"TypeB" (40 cells). Names correspond to individual cell barcodes
(`"Cell_01"` to `"Cell_80"`).

## Source

Derived from `example_scrna$cell_types`.

## Examples

``` r
data(example_cell_types)
head(example_cell_types)
#> Cell_01 Cell_02 Cell_03 Cell_04 Cell_05 Cell_06 
#>   TypeA   TypeA   TypeA   TypeA   TypeA   TypeA 
#> Levels: TypeA TypeB
table(example_cell_types)
#> example_cell_types
#> TypeA TypeB 
#>    40    40 
```
