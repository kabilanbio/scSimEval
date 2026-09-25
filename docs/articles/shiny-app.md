# Interactive Web Application (Shiny Studio)

## Introduction

For researchers who prefer an interactive graphical interface rather
than writing R code, **`scSimEval`** includes an easy-to-use, embedded
Shiny web application:
**[`launch_scSimEval_app()`](https://kabilanbio.github.io/scSimEval/reference/launch_scSimEval_app.md)**.

The Shiny application allows you to explore benchmark results across all
62 evaluation measures, adjust category weights to re-rank simulation
methods, and export high-resolution (600 DPI) figures directly in your
web browser.

------------------------------------------------------------------------

## Launching the Application

To start the interactive studio, load `scSimEval` and run:

``` r
library(scSimEval)

# Launch the app locally in your default web browser
launch_scSimEval_app()
```

By default, the application opens an interactive session in your browser
(at `http://127.0.0.1:port`).

------------------------------------------------------------------------

## Key Features

#### 1. Instant Demo Benchmark Mode

- Includes pre-computed benchmark results across multiple single-cell
  simulators with zero waiting time.
- Allows immediate exploration of bubble matrices, category rankings,
  and diagnostic plots without needing to execute simulations first.

#### 2. Custom Benchmark Upload

- Upload custom evaluation outputs generated from
  [`evaluate_multiomics_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_accuracy.md),
  [`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md),
  or
  [`evaluate_multiple_datasets()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiple_datasets.md).
- Automatically checks table structures and organizes metric categories.

#### 3. Dynamic Category Weighting & Leaderboard Re-Ranking

- Customize the relative importance of each of the 8 categories using
  interactive sliders.
- Observe real-time changes to composite simulator scores and
  leaderboard rankings based on your specific priorities
  (e.g. prioritizing marker genes over computer speed).

#### 4. Interactive Bubble Matrix Explorer

- Filter metrics by category or specific biological properties.
- Inspect exact standardized scores ($`0.00`$ to $`1.00`$) and
  unnormalized raw values on hover.
- Adjust bubble scaling and highlight threshold glyphs interactively.

#### 5. High-Resolution (600 DPI) Figure Export

- Export benchmark figures directly from the interface in
  high-resolution formats (PNG, PDF, SVG) at **600 DPI**.
- Download consolidated benchmark summary tables in CSV and Excel
  formats.

------------------------------------------------------------------------

## Summary

The interactive Shiny Studio provides an accessible, zero-code
environment to explore, customize, and communicate single-cell
simulation benchmarking results.
