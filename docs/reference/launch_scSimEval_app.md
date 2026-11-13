# Launch Interactive scSimEval Benchmarking Studio

Launches the Shiny web application embedded in the scSimEval package.
Provides an intuitive graphical interface to upload biological reference
and simulated datasets (single-cell scRNA-seq, scATAC-seq, or paired
multiomics), specify computational scalability metrics (elapsed time and
peak memory), interactively inspect the flagship 48-measure comparative
bubble matrix, explore diagnostic figures, and export high-resolution
figures (600 DPI publication quality), multi-page PDF reports, Excel
(.xlsx) workbooks, and complete results ZIP packages.

## Usage

``` r
launch_scSimEval_app(
  port = NULL,
  host = "127.0.0.1",
  launch.browser = interactive(),
  maxRequestSize = 10 * 1024^3
)
```

## Arguments

- port:

  Optional port number for the local web server. Default is `NULL`
  (random open port).

- host:

  Character string specifying the IP address to listen on. Defaults to
  `"127.0.0.1"`.

- launch.browser:

  Logical, whether to automatically launch the default web browser.
  Defaults to `TRUE` in interactive sessions.

- maxRequestSize:

  Numeric value specifying the maximum allowable upload file size in
  bytes. Defaults to `10 * 1024^3` (10 GB) to comfortably accommodate
  large single-cell multiomics matrices.

## Value

Invisibly returns the Shiny app process object.

## Examples

``` r
# Locate the embedded Shiny app directory bundled with the package
system.file("shiny", "scSimEvalApp", package = "scSimEval")
#> [1] "C:/Users/kabil/AppData/Local/R/win-library/4.6/scSimEval/shiny/scSimEvalApp"
# \donttest{
# Launch the interactive studio (interactive sessions only)
if (interactive()) {
  launch_scSimEval_app()
}
# }
```
