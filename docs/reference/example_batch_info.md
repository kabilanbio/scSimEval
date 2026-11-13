# Example Batch Annotations Vector

A named factor vector containing technical batch annotations ("Batch1",
"Batch2") for 80 single cells, matching the column names of
[`example_scrna`](https://kabilanbio.github.io/scSimEval/reference/example_scrna.md)
and
[`example_scatac`](https://kabilanbio.github.io/scSimEval/reference/example_scatac.md).
Useful as a reference template for formatting batch metadata inputs in
`scSimEval`.

## Usage

``` r
data(example_batch_info)
```

## Format

A named factor with 80 elements across 2 levels: "Batch1" (40 cells) and
"Batch2" (40 cells). Names correspond to individual cell barcodes
(`"Cell_01"` to `"Cell_80"`).

## Source

Derived from `example_scrna$batch_info`.

## Examples

``` r
data(example_batch_info)
head(example_batch_info)
#> Cell_01 Cell_02 Cell_03 Cell_04 Cell_05 Cell_06 
#>  Batch1  Batch2  Batch1  Batch2  Batch1  Batch2 
#> Levels: Batch1 Batch2
table(example_batch_info)
#> example_batch_info
#> Batch1 Batch2 
#>     40     40 
```
