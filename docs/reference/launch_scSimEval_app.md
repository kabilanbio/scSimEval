# Launch Interactive scSimEval Benchmarking Studio

Launches the Shiny web application embedded in the scSimEval package.
Provides an intuitive graphical interface to ingest single-cell
reference and simulated datasets, configure ground-truth-free evaluation
pipelines, interactively inspect the flagship 62-measure bubble matrix,
adjust category weighting for custom method rankings, and export
high-resolution (600 DPI) figures and tables.

## Usage

``` r
launch_scSimEval_app(
  port = NULL,
  host = "127.0.0.1",
  launch.browser = interactive()
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

## Value

Invisibly returns the Shiny app process object.

## Examples

``` r
if (FALSE) { # \dontrun{
library(scSimEval)
# Launch the interactive studio
launch_scSimEval_app()

} # }
```
