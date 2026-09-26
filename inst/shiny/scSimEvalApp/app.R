# inst/shiny/scSimEvalApp/app.R
# Unified Single-Cell & Multiomics Simulation Benchmarking Studio
# Powered by scSimEval (62 Curated Ground-Truth-Free Measures)

library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(Matrix)
library(scSimEval)

# Set max upload size to 500 MB for large single-cell datasets
options(shiny.maxRequestSize = 500 * 1024^2)

# Load demo benchmark data if present
demo_data_path <- system.file("shiny", "scSimEvalApp", "data", "demo_benchmark_data.rds", package = "scSimEval")
if (demo_data_path == "" || !file.exists(demo_data_path)) {
  demo_data_path <- file.path("data", "demo_benchmark_data.rds")
}
initial_demo <- if (file.exists(demo_data_path)) readRDS(demo_data_path) else NULL

# ==============================================================================
# Helper Functions: Robust Matrix and Label Ingestion
# ==============================================================================
read_uploaded_matrix <- function(file_path, file_name) {
  if (is.null(file_path) || !file.exists(file_path)) return(NULL)
  ext <- tolower(tools::file_ext(file_name))
  
  if (ext == "rds") {
    obj <- readRDS(file_path)
    if (is.matrix(obj) || inherits(obj, "Matrix") || is.data.frame(obj)) {
      return(as.matrix(obj))
    } else if (inherits(obj, "SingleCellExperiment") && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      return(as.matrix(SingleCellExperiment::counts(obj)))
    } else if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
      return(as.matrix(Seurat::GetAssayData(obj, slot = "counts")))
    } else if (is.list(obj) && length(obj) > 0 && (is.matrix(obj[[1]]) || is.data.frame(obj[[1]]))) {
      return(as.matrix(obj[[1]]))
    } else {
      stop("Unsupported RDS file structure. Please upload a matrix or data.frame.")
    }
  } else if (ext == "csv") {
    df <- utils::read.csv(file_path, row.names = 1, check.names = FALSE)
    return(as.matrix(df))
  } else if (ext %in% c("tsv", "txt")) {
    df <- utils::read.table(file_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    return(as.matrix(df))
  } else {
    stop(paste("Unsupported file format:", ext))
  }
}

read_uploaded_labels <- function(file_path, file_name) {
  if (is.null(file_path) || !file.exists(file_path)) return(NULL)
  ext <- tolower(tools::file_ext(file_name))
  if (ext == "rds") {
    obj <- readRDS(file_path)
    if (is.vector(obj) || is.factor(obj)) return(as.factor(obj))
    if (is.data.frame(obj)) return(as.factor(obj[[1]]))
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- if (ext == "csv") "," else "\t"
    df <- utils::read.table(file_path, sep = sep, header = TRUE, stringsAsFactors = FALSE)
    return(as.factor(df[[1]]))
  }
  return(NULL)
}

# Format Excel workbook
export_excel_workbook <- function(file, benchmark_df, leaderboard_df = NULL) {
  sheets <- list(All_Benchmark_Metrics = benchmark_df)
  if (!is.null(leaderboard_df) && nrow(leaderboard_df) > 0) {
    sheets$Method_Rankings <- leaderboard_df
  }
  sc_df <- subset(benchmark_df, grepl("Scalability", Category))
  if (nrow(sc_df) > 0) {
    sheets$Scalability_Metrics <- sc_df
  }
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(sheets, path = file)
  } else if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(sheets, file = file)
  } else {
    utils::write.csv(benchmark_df, file, row.names = FALSE)
  }
}

# Export High-Res JPEG
export_single_jpeg <- function(file, plot_obj, width = 14, height = 9, dpi = 300) {
  ggplot2::ggsave(file, plot = plot_obj, device = "jpeg", width = width, height = height, dpi = dpi)
}

# Generate Multi-Page PDF Report
generate_all_plots_pdf <- function(file, benchmark_data, toy_ref = NULL, toy_sim = NULL) {
  grDevices::pdf(file, width = 14, height = 9, onefile = TRUE)
  try(print(plot_benchmark_bubble_matrix(benchmark_data, base_size = 9)), silent = TRUE)
  try(print(plot_evaluation_summary(benchmark_data, base_size = 10)), silent = TRUE)
  try(print(plot_scalability_benchmark(benchmark_data, base_size = 10)), silent = TRUE)
  try(print(plot_metric_boxplots(benchmark_data, base_size = 10)), silent = TRUE)
  try(print(plot_metric_heatmap(benchmark_data, base_size = 10)), silent = TRUE)
  try(print(plot_metric_pca(benchmark_data, base_size = 10)), silent = TRUE)
  try(print(plot_metric_mds(benchmark_data, base_size = 10)), silent = TRUE)
  if (!is.null(toy_ref) && !is.null(toy_sim)) {
    try(print(plot_distribution_qc(toy_ref, toy_sim, base_size = 9)), silent = TRUE)
  }
  grDevices::dev.off()
}

# ==============================================================================
# Custom Theme & Aesthetics
# ==============================================================================
app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#1B4F72",
  secondary = "#16A085",
  success = "#27AE60",
  info = "#2980B9",
  warning = "#E67E22",
  danger = "#C0392B",
  base_font = font_google("Inter")
)

# ==============================================================================
# UI Definition
# ==============================================================================
ui <- page_navbar(
  title = "scSimEval Studio",
  id = "nav_active",
  theme = app_theme,
  fillable = TRUE,
  
  header = tags$head(
    tags$style(HTML("
      .navbar { box-shadow: 0 2px 8px rgba(0,0,0,0.1); font-weight: 600; }
      .nav-link { font-size: 0.95rem; }
      .stat-card { border-radius: 10px; border-left: 5px solid #1B4F72; box-shadow: 0 3px 6px rgba(0,0,0,0.05); background: white; padding: 16px; margin-bottom: 15px; }
      .stat-number { font-size: 2.2rem; font-weight: 800; color: #1B4F72; line-height: 1; }
      .stat-label { font-size: 0.82rem; text-transform: uppercase; color: #7F8C8D; font-weight: 600; letter-spacing: 0.5px; }
      .category-pill { display: inline-block; padding: 5px 12px; border-radius: 14px; font-size: 0.82rem; font-weight: 600; color: white; margin: 3px; }
      .hero-box { background: linear-gradient(135deg, #1B4F72 0%, #2C3E50 100%); color: white; border-radius: 12px; padding: 26px; margin-bottom: 20px; box-shadow: 0 4px 12px rgba(27,79,114,0.25); }
      .card-header { font-weight: 700; color: #2C3E50; background-color: #F8F9F9; }
      .btn-primary { background-color: #1B4F72; border-color: #1B4F72; }
      .btn-primary:hover { background-color: #154360; border-color: #154360; }
      .btn-success { background-color: #16A085; border-color: #16A085; }
      .btn-success:hover { background-color: #117A65; border-color: #117A65; }
      .guide-step { background: #FFFFFF; border-radius: 10px; border: 1px solid #E2E8F0; padding: 18px; margin-bottom: 15px; }
      .guide-num { display: inline-block; width: 32px; height: 32px; line-height: 32px; border-radius: 50%; background: #1B4F72; color: white; font-weight: 800; text-align: center; margin-right: 10px; }
      .bubble-plot-box { background: white; border-radius: 8px; padding: 12px; box-shadow: 0 2px 6px rgba(0,0,0,0.06); overflow-x: auto; }
    "))
  ),
  
  # ============================================================================
  # TAB 1: OVERVIEW & USER GUIDE
  # ============================================================================
  nav_panel(
    "Overview & Guide",
    fluidRow(
      column(
        12,
        div(
          class = "hero-box",
          h2("scSimEval: Single-Cell & Multiomics Simulation Benchmarking Studio", style = "font-weight: 800;"),
          p("A simple, unified evaluation toolkit to benchmark how realistically simulated data mirror real biological datasets. Evaluates single-cell and multiomics simulation methods across 62 curated metrics and 8 canonical categories without requiring hardcoded synthetic ground truth.", style = "font-size: 1.05rem; opacity: 0.95;"),
          hr(style = "border-color: rgba(255,255,255,0.25);"),
          div(
            actionButton("btn_go_data", "Step 1: Load or Ingest Data", class = "btn btn-outline-light me-2", icon = icon("upload")),
            actionButton("btn_go_bubble", "Step 2: Explore Comparative Bubble Matrix", class = "btn btn-success me-2", icon = icon("chart-pie")),
            actionButton("btn_go_viz", "Step 3: View All 8 Diagnostic Figures", class = "btn btn-info me-2 text-white", icon = icon("images")),
            actionButton("btn_go_download", "Step 4: Download Reports & Excel (.xlsx)", class = "btn btn-outline-light", icon = icon("download"))
          )
        )
      )
    ),
    fluidRow(
      column(3, div(class = "stat-card", div(class = "stat-number", "62"), div(class = "stat-label", "Curated Evaluation Metrics"))),
      column(3, div(class = "stat-card", div(class = "stat-number", "8"), div(class = "stat-label", "Canonical Biological Categories"))),
      column(3, div(class = "stat-card", div(class = "stat-number", "2"), div(class = "stat-label", "Scalability Metrics (Runtime & RAM)"))),
      column(3, div(class = "stat-card", div(class = "stat-number", "100%"), div(class = "stat-label", "Ground-Truth-Free Evaluation")))
    ),
    fluidRow(
      column(
        7,
        card(
          card_header("Eight Canonical Evaluation Categories"),
          card_body(
            tags$div(
              tags$span(class = "category-pill", style = "background-color: #2E86AB;", "(I) Distributional Properties (14 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #17A589;", "(II) Correlations & Zero-Inflation (6 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #C0392B;", "(III) Cellular Structure & Concordance (10 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #D4AC0D;", "(IV) Batch Effects & Confounder Mixing (7 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #CA6F1E;", "(V) Biological Signal & Downstream Fidelity (7 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #7D3C98;", "(VI) Trajectory & Lineage Dynamics (2 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #2E4057;", "(VII) Cross-Modal Coupling & Modularity (6 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #1E8449;", "(VIII) Computational Scalability (2 metrics)")
            ),
            hr(),
            h5("Simple, Ground-Truth-Free Methodology", style = "font-weight: 700;"),
            p("Standard benchmarking often depends on pre-defined synthetic labels (such as hardcoded differentially expressed genes or preset gene regulatory networks). This can cause circular validation where simulators are scored on their own built-in assumptions."),
            p("In contrast, ", tags$b("scSimEval"), " directly measures statistical differences, cellular clustering consistency, and manifold geometry against real biological reference datasets.")
          )
        )
      ),
      column(
        5,
        card(
          card_header("Beginner's Step-by-Step Workflow"),
          card_body(
            div(class = "guide-step",
                div(class = "guide-num", "1"),
                tags$b("Ingest Data: "),
                "Click 'Load 6-Simulator Demo Benchmark' to immediately explore pre-calculated results, or upload your own real reference and simulated datasets (individual or multiple simulators)."
            ),
            div(class = "guide-step",
                div(class = "guide-num", "2"),
                tags$b("Comparative Bubble Matrix: "),
                "Compare all simulators side-by-side in one comprehensive matrix. Larger bubbles indicate superior fidelity to the biological reference."
            ),
            div(class = "guide-step",
                div(class = "guide-num", "3"),
                tags$b("Download Complete Results: "),
                "Export results as Excel (.xlsx), CSV, RDS, high-resolution JPEGs (300 DPI), multi-page PDF report, or download everything as a single .zip file."
            )
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 2: DATA HUB (UPLOAD & EVALUATE)
  # ============================================================================
  nav_panel(
    "Data Hub",
    layout_sidebar(
      sidebar = sidebar(
        width = 380,
        title = "Dataset Ingestion & Live Evaluation",
        
        # Ingestion Mode Selector
        radioButtons(
          "opt_data_mode", "Choose Ingestion Mode:",
          choices = c(
            "Mode 1: Explore Demo Benchmark (6 Simulators)" = "demo",
            "Mode 2: Evaluate Single Simulator" = "single",
            "Mode 3: Evaluate Multiple Simulators (Batch)" = "multi",
            "Mode 4: Upload Pre-Computed Benchmark (.rds / .csv)" = "upload_bench"
          ),
          selected = "demo"
        ),
        hr(),
        
        # ----------------- Mode 1: Demo -----------------
        conditionalPanel(
          condition = "input.opt_data_mode == 'demo'",
          p("Instantly explore 6 leading simulators (Splatter, scDesign3, SCRIP, SymSim, dyngen, simATAC) across all 62 evaluation measures.", style = "font-size: 0.9rem; color: #555;"),
          actionButton("btn_load_demo", "Load 6-Simulator Demo Benchmark", class = "btn btn-success w-100", icon = icon("play"))
        ),
        
        # ----------------- Mode 2: Single Simulator -----------------
        conditionalPanel(
          condition = "input.opt_data_mode == 'single'",
          h6(tags$b("1. Upload Count Matrices")),
          fileInput("file_single_ref", "Reference Count Matrix (.rds / .csv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
          fileInput("file_single_sim", "Simulated Count Matrix (.rds / .csv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
          textInput("txt_single_name", "Simulator Method Name:", value = "MySimulator"),
          
          h6(tags$b("2. Computational Scalability (2 Metrics)")),
          fluidRow(
            column(6, numericInput("num_single_time", "Elapsed Time (s):", value = 35.0, min = 0.1, step = 0.5)),
            column(6, numericInput("num_single_mem", "Peak RAM (MB):", value = 820.0, min = 1, step = 10))
          ),
          
          h6(tags$b("3. Optional Annotations")),
          fileInput("file_single_celltypes", "Optional Cell Types (.rds / .csv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
          fileInput("file_single_batch", "Optional Batch Labels (.rds / .csv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
          
          checkboxInput("chk_append_single", "Append to existing benchmark (compare together)", value = TRUE),
          actionButton("btn_run_single_eval", "Run Single Evaluation", class = "btn btn-primary w-100", icon = icon("calculator"))
        ),
        
        # ----------------- Mode 3: Multiple Simulators -----------------
        conditionalPanel(
          condition = "input.opt_data_mode == 'multi'",
          h6(tags$b("1. Reference Biological Dataset")),
          fileInput("file_multi_ref", "Reference Count Matrix (.rds / .csv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
          
          h6(tags$b("2. Multiple Simulated Datasets")),
          fileInput("file_multi_sims", "Simulated Count Matrices (Select 1 or more):", multiple = TRUE, accept = c(".rds", ".csv", ".tsv", ".txt")),
          p("After selecting simulated files, configure their simulator names and scalability metrics below:", style = "font-size: 0.85rem; color: #555;"),
          uiOutput("ui_multi_sim_scalability_inputs"),
          
          hr(),
          actionButton("btn_run_multi_eval", "Run Multi-Simulator Evaluation", class = "btn btn-primary w-100", icon = icon("cogs"))
        ),
        
        # ----------------- Mode 4: Pre-Computed Benchmark -----------------
        conditionalPanel(
          condition = "input.opt_data_mode == 'upload_bench'",
          p("Upload benchmark results previously saved from scSimEval evaluation pipelines.", style = "font-size: 0.9rem; color: #555;"),
          fileInput("file_bench_upload", "Upload scSimEval Object (.rds or .csv):", accept = c(".rds", ".csv")),
          actionButton("btn_load_uploaded_bench", "Load Uploaded Benchmark", class = "btn btn-info text-white w-100", icon = icon("folder-open"))
        )
      ),
      
      # Main Card: Status & Data Table Preview
      card(
        card_header("Active Benchmark Dataset Status"),
        card_body(
          uiOutput("ui_status_banner"),
          hr(),
          h5("Benchmark Summary Table Preview"),
          p("Showing evaluated metrics across simulators. Column 'Score' represents direction-aware standardized fidelity in [0, 1].", style = "font-size: 0.88rem; color: #7F8C8D;"),
          DTOutput("table_active_data_preview")
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 3: COMPARATIVE BUBBLE MATRIX (FLAGSHIP FIGURE)
  # ============================================================================
  nav_panel(
    "Comparative Bubble Matrix",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        title = "Display & Ranking Controls",
        
        selectInput(
          "sel_bubble_cat", "Filter by Category:",
          choices = c(
            "All Categories (All 62 Measures)" = "all",
            "(I) Distributional Properties",
            "(II) Correlations & Zero-Inflation",
            "(III) Cellular Structure & Concordance",
            "(IV) Batch Effects & Confounder Mixing",
            "(V) Biological Signal & Downstream Fidelity",
            "(VI) Trajectory & Lineage Dynamics",
            "(VII) Cross-Modal Coupling & Modularity",
            "(VIII) Computational Scalability"
          ),
          selected = "all"
        ),
        
        uiOutput("ui_bubble_method_picker"),
        
        sliderInput("sld_bubble_size", "Bubble Size Range:", min = 1, max = 16, value = c(2.5, 8.5), step = 0.5),
        checkboxInput("chk_bubble_missing", "Show Missing Dots (Unmeasured)", value = TRUE),
        checkboxInput("chk_bubble_norm", "Normalize Scores to [0, 1]", value = TRUE),
        
        selectInput(
          "sel_bubble_height", "Plot Display Height (Bigger Size):",
          choices = c(
            "Standard (850px)" = "850px",
            "Large (1000px)" = "1000px",
            "Extra Large (1250px)" = "1250px"
          ),
          selected = "850px"
        ),
        
        hr(),
        h6(tags$b("Custom Category Weights for Leaderboard:")),
        sliderInput("wt_cat1", "Cat I (Distribution):", min = 0, max = 3, value = 1, step = 0.5),
        sliderInput("wt_cat3", "Cat III (Cellular Structure):", min = 0, max = 3, value = 1, step = 0.5),
        sliderInput("wt_cat5", "Cat V (Biological Signals):", min = 0, max = 3, value = 1, step = 0.5),
        sliderInput("wt_cat8", "Cat VIII (Scalability):", min = 0, max = 3, value = 1, step = 0.5),
        
        hr(),
        h6(tags$b("Download This Plot:")),
        downloadButton("download_bubble_jpeg", "Download JPEG (300 DPI)", class = "btn btn-primary w-100 mb-2"),
        downloadButton("download_bubble_pdf", "Download Vector PDF", class = "btn btn-outline-secondary w-100")
      ),
      
      card(
        card_header("Multi-Dimensional Simulation Fidelity Bubble Matrix (All Datasets Compared)"),
        card_body(
          p("Each column represents a simulator method; each row represents a curated evaluation metric. Bubble size reflects standardized fidelity (larger bubbles = higher fidelity to reference). Color indicates biological category.", style = "font-size: 0.9rem; color: #555;"),
          div(class = "bubble-plot-box",
              uiOutput("ui_bubble_plot_container")
          ),
          hr(),
          h5("Method Ranking Leaderboard (Weighted Overall Fidelity Score)", style = "font-weight: 700;"),
          DTOutput("table_bubble_leaderboard")
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 4: DIAGNOSTIC VISUALIZATIONS (ALL 8 FIGURES)
  # ============================================================================
  nav_panel(
    "Visualizations",
    navset_pill(
      # Sub-tab 1: Overall Evaluation Summary
      nav_panel(
        "1. Evaluation Summary",
        card(
          card_header("Overall Category Evaluation Summary (plot_evaluation_summary)"),
          card_body(
            fluidRow(
              column(4, checkboxInput("chk_sum_labels", "Show Score Labels", value = TRUE)),
              column(4, checkboxInput("chk_sum_norm", "Normalize Scores [0, 1]", value = TRUE)),
              column(4,
                     downloadButton("download_sum_jpeg", "Download JPEG", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_sum_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_eval_summary", height = "520px")
          )
        )
      ),
      
      # Sub-tab 2: Distribution QC
      nav_panel(
        "2. Distribution QC",
        card(
          card_header("Empirical vs Simulated Distribution Quality (plot_distribution_qc)"),
          card_body(
            fluidRow(
              column(4, selectInput("sel_dist_layout", "QC Layout:", choices = c("Comprehensive" = "comprehensive", "Density Curves" = "density"), selected = "comprehensive")),
              column(4, p("Compares expression densities, library sizes, and zero-inflation between real reference and simulated cells.", style = "font-size: 0.85rem; color: #666;")),
              column(4,
                     downloadButton("download_dist_jpeg", "Download JPEG", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_dist_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_dist_qc", height = "550px")
          )
        )
      ),
      
      # Sub-tab 3: Scalability Benchmark
      nav_panel(
        "3. Scalability Benchmark",
        card(
          card_header("Computational Scalability Benchmark (plot_scalability_benchmark)"),
          card_body(
            fluidRow(
              column(4,
                     selectInput(
                       "sel_scale_type", "Scalability View:",
                       choices = c(
                         "4-Panel Comprehensive Layout" = "composite",
                         "Runtime Execution Time" = "runtime",
                         "Peak RAM Memory Usage" = "memory",
                         "Runtime vs Memory Trade-Off" = "tradeoff",
                         "Resource Cost Footprint" = "cost",
                         "Cell Throughput" = "throughput"
                       ),
                       selected = "composite"
                     )
              ),
              column(4, sliderInput("sld_scale_cells", "Hypothetical Cell Count (for Throughput):", min = 500, max = 50000, value = 5000, step = 500)),
              column(4,
                     downloadButton("download_scale_jpeg", "Download JPEG", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_scale_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_scale_bench", height = "600px")
          )
        )
      ),
      
      # Sub-tab 4: Metric Boxplots
      nav_panel(
        "4. Metric Boxplots",
        card(
          card_header("Fidelity Score Distributions by Category (plot_metric_boxplots)"),
          card_body(
            fluidRow(
              column(4, selectInput("sel_box_score_type", "Score Type:", choices = c("Normalized [0, 1]" = "normalized", "Raw Metric Value" = "raw"), selected = "normalized")),
              column(4, selectInput("sel_box_facet", "Facet By:", choices = c("Category" = "category", "Metric" = "metric"), selected = "category")),
              column(4,
                     downloadButton("download_box_jpeg", "Download JPEG", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_box_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_metric_boxes", height = "550px")
          )
        )
      ),
      
      # Sub-tab 5: Metric Heatmap
      nav_panel(
        "5. Metric Heatmap",
        card(
          card_header("Cross-Metric Correlation & Performance Heatmap (plot_metric_heatmap)"),
          card_body(
            fluidRow(
              column(4, checkboxInput("chk_heat_cluster_rows", "Cluster Simulators (Rows)", value = FALSE)),
              column(4, checkboxInput("chk_heat_cluster_cols", "Cluster Metrics (Columns)", value = FALSE)),
              column(4,
                     downloadButton("download_heat_jpeg", "Download JPEG", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_heat_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_metric_heat", height = "550px")
          )
        )
      ),
      
      # Sub-tab 6: PCA Ordination
      nav_panel(
        "6. PCA Ordination",
        card(
          card_header("PCA Projection of Simulators in Metric Performance Space (plot_metric_pca)"),
          card_body(
            fluidRow(
              column(4, selectInput("sel_pca_panel", "Panel View:", choices = c("Both (Biplot + Loadings)" = "both", "Simulators Only" = "methods", "Loadings Only" = "loadings"), selected = "both")),
              column(4, sliderInput("sld_pca_loadings", "Top Metric Loadings to Display:", min = 4, max = 20, value = 10, step = 1)),
              column(4,
                     downloadButton("download_pca_jpeg", "Download JPEG", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_pca_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_metric_pca_out", height = "520px")
          )
        )
      ),
      
      # Sub-tab 7: MDS Ordination
      nav_panel(
        "7. MDS Metric Space",
        card(
          card_header("Multi-Dimensional Scaling (MDS) Ordination (plot_metric_mds)"),
          card_body(
            fluidRow(
              column(4, selectInput("sel_mds_by", "MDS Ordination Target:", choices = c("By Simulators" = "simulators", "By Metric Summaries" = "summaries"), selected = "simulators")),
              column(4, p("Maps simulators onto a 2D Euclidean distance space derived from multi-category fidelity scores.", style = "font-size: 0.85rem; color: #666;")),
              column(4,
                     downloadButton("download_mds_jpeg", "Download JPEG", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_mds_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_metric_mds_out", height = "520px")
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 5: DOWNLOAD & REPORT CENTER
  # ============================================================================
  nav_panel(
    "Download Center",
    fluidRow(
      column(
        4,
        card(
          card_header("1. Complete Results Package (.zip)"),
          card_body(
            p("Download an all-in-one ZIP archive containing:"),
            tags$ul(
              tags$li(tags$b("Excel Workbook (.xlsx): "), "All metrics, method leaderboard, and scalability tables."),
              tags$li(tags$b("Tidy CSV Table (.csv): "), "Full 62-metric evaluation matrix."),
              tags$li(tags$b("R Data Object (.rds): "), "Reproducible scSimEval benchmark object."),
              tags$li(tags$b("Multi-Page PDF Report: "), "High-resolution PDF containing all 8 active figures."),
              tags$li(tags$b("Figure Images: "), "Individual 300 DPI publication-ready JPEGs.")
            ),
            hr(),
            downloadButton("download_complete_zip", "Download Complete Results (.zip)", class = "btn btn-success w-100 py-2")
          )
        )
      ),
      column(
        4,
        card(
          card_header("2. Tabular Data Formats"),
          card_body(
            p("Export the active benchmark data table for supplementary materials or statistical analysis:"),
            downloadButton("download_excel", "Download All Metrics (.xlsx)", class = "btn btn-primary w-100 mb-2"),
            downloadButton("download_csv", "Download Master Table (.csv)", class = "btn btn-outline-primary w-100 mb-2"),
            downloadButton("download_rds", "Download Results Object (.rds)", class = "btn btn-outline-secondary w-100")
          )
        )
      ),
      column(
        4,
        card(
          card_header("3. Multi-Page PDF Figure Report"),
          card_body(
            p("Generate a comprehensive multi-page PDF compilation containing all 8 evaluation figures:"),
            tags$ol(
              tags$li("Flagship 62-Metric Bubble Matrix"),
              tags$li("Category Evaluation Summary"),
              tags$li("Computational Scalability Benchmark"),
              tags$li("Metric Boxplots & Variance"),
              tags$li("Metric Correlation Heatmap"),
              tags$li("PCA Simulator Ordination"),
              tags$li("MDS Metric Space"),
              tags$li("Distribution QC Curves")
            ),
            hr(),
            downloadButton("download_all_plots_pdf", "Download All Figures (.pdf)", class = "btn btn-info text-white w-100 py-2")
          )
        )
      )
    ),
    fluidRow(
      column(
        12,
        card(
          card_header("Interactive Master Benchmark Data Table"),
          card_body(
            DTOutput("table_master_export")
          )
        )
      )
    )
  )
)

# ==============================================================================
# SERVER DEFINITION
# ==============================================================================
server <- function(input, output, session) {
  
  # Reactive state container
  rv <- reactiveValues(
    benchmark_df = if (!is.null(initial_demo)) initial_demo$benchmark_summary_table else NULL,
    methods = if (!is.null(initial_demo)) initial_demo$methods else NULL,
    toy_ref = if (!is.null(initial_demo)) initial_demo$toy_data$ref else NULL,
    toy_sim = if (!is.null(initial_demo)) initial_demo$toy_data$sim else NULL,
    source_name = if (!is.null(initial_demo)) "Built-in 6-Simulator Demo Benchmark" else "No Data Loaded"
  )
  
  # Quick navigation triggers
  observeEvent(input$btn_go_data, { nav_select("nav_active", "Data Hub") })
  observeEvent(input$btn_go_bubble, { nav_select("nav_active", "Comparative Bubble Matrix") })
  observeEvent(input$btn_go_viz, { nav_select("nav_active", "Visualizations") })
  observeEvent(input$btn_go_download, { nav_select("nav_active", "Download Center") })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 1 - Load Demo Benchmark
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_load_demo, {
    if (!is.null(initial_demo)) {
      rv$benchmark_df <- initial_demo$benchmark_summary_table
      rv$methods <- initial_demo$methods
      rv$toy_ref <- initial_demo$toy_data$ref
      rv$toy_sim <- initial_demo$toy_data$sim
      rv$source_name <- "Built-in 6-Simulator Demo Benchmark (Splatter, scDesign3, SCRIP, SymSim, dyngen, simATAC)"
      
      updateCheckboxGroupInput(session, "sel_bubble_methods", choices = rv$methods, selected = rv$methods)
      showNotification("Successfully loaded 6-Simulator Demo Benchmark!", type = "message")
    } else {
      showNotification("Demo benchmark file not found on disk.", type = "warning")
    }
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Dynamic Scalability Inputs for Multiple Simulators Upload (Mode 3)
  # ----------------------------------------------------------------------------
  output$ui_multi_sim_scalability_inputs <- renderUI({
    req(input$file_multi_sims)
    n_files <- nrow(input$file_multi_sims)
    if (n_files == 0) return(NULL)
    
    inputs_list <- lapply(seq_len(n_files), function(i) {
      fname <- input$file_multi_sims$name[i]
      default_name <- tools::file_path_sans_ext(fname)
      
      div(
        style = "background: #F8F9FA; border-left: 3px solid #1B4F72; padding: 10px; margin-bottom: 8px; border-radius: 4px;",
        tags$b(paste0("Simulator ", i, ": "), style = "font-size: 0.9rem; color: #1B4F72;"),
        textInput(paste0("multi_sim_name_", i), "Method Name:", value = default_name),
        fluidRow(
          column(6, numericInput(paste0("multi_sim_time_", i), "Elapsed Time (s):", value = round(25 + i * 12, 1), min = 0.1, step = 0.5)),
          column(6, numericInput(paste0("multi_sim_mem_", i), "Peak RAM (MB):", value = round(650 + i * 180, 0), min = 1, step = 10))
        )
      )
    })
    
    tagList(inputs_list)
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 2 - Run Single Simulator Evaluation
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_run_single_eval, {
    req(input$file_single_ref, input$file_single_sim)
    
    withProgress(message = "Evaluating simulator...", detail = "Reading count matrices", value = 0.2, {
      tryCatch({
        ref_mat <- read_uploaded_matrix(input$file_single_ref$datapath, input$file_single_ref$name)
        sim_mat <- read_uploaded_matrix(input$file_single_sim$datapath, input$file_single_sim$name)
        
        incProgress(0.3, detail = "Computing 62 evaluation metrics...")
        
        sim_name <- trimws(input$txt_single_name)
        if (sim_name == "") sim_name <- "Simulator"
        
        elapsed_sec <- as.numeric(input$num_single_time)
        peak_ram <- as.numeric(input$num_single_mem)
        
        res <- evaluate_simulation_accuracy(
          ref_data = ref_mat,
          sim_data = sim_mat,
          memory_mb = peak_ram,
          elapsed_time = elapsed_sec,
          compute_bivariate = FALSE,
          verbose = FALSE
        )
        
        incProgress(0.3, detail = "Formatting results table...")
        
        tbl <- res$metrics_summary_table
        tbl$Method <- sim_name
        
        # Harmonize column names
        if (!"Score" %in% colnames(tbl) && "Value" %in% colnames(tbl)) {
          tbl$Score <- tbl$Value
        }
        
        if (isTRUE(input$chk_append_single) && !is.null(rv$benchmark_df)) {
          # Remove any existing rows for this simulator name
          existing_clean <- rv$benchmark_df[rv$benchmark_df$Method != sim_name, , drop = FALSE]
          combined <- rbind(existing_clean, tbl[, intersect(colnames(existing_clean), colnames(tbl))])
          rv$benchmark_df <- combined
          rv$methods <- unique(combined$Method)
          rv$source_name <- paste0("Combined Benchmark (", length(rv$methods), " Simulators)")
        } else {
          rv$benchmark_df <- tbl
          rv$methods <- sim_name
          rv$source_name <- paste0("Single Evaluation: ", sim_name)
        }
        
        rv$toy_ref <- ref_mat
        rv$toy_sim <- sim_mat
        
        updateCheckboxGroupInput(session, "sel_bubble_methods", choices = rv$methods, selected = rv$methods)
        incProgress(0.2, detail = "Done!")
        showNotification(paste0("Successfully evaluated '", sim_name, "'!"), type = "message")
      }, error = function(e) {
        showNotification(paste("Evaluation error:", e$message), type = "error")
      })
    })
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 3 - Run Multiple Simulators Evaluation
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_run_multi_eval, {
    req(input$file_multi_ref, input$file_multi_sims)
    n_files <- nrow(input$file_multi_sims)
    req(n_files > 0)
    
    withProgress(message = "Batch Simulators Evaluation", value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Loading biological reference matrix...")
        ref_mat <- read_uploaded_matrix(input$file_multi_ref$datapath, input$file_multi_ref$name)
        
        results_list <- list()
        first_sim_mat <- NULL
        
        for (i in seq_len(n_files)) {
          sim_name <- input[[paste0("multi_sim_name_", i)]]
          if (is.null(sim_name) || trimws(sim_name) == "") {
            sim_name <- tools::file_path_sans_ext(input$file_multi_sims$name[i])
          }
          
          sim_time <- as.numeric(input[[paste0("multi_sim_time_", i)]])
          if (is.null(sim_time) || is.na(sim_time)) sim_time <- 30.0
          
          sim_mem <- as.numeric(input[[paste0("multi_sim_mem_", i)]])
          if (is.null(sim_mem) || is.na(sim_mem)) sim_mem <- 800.0
          
          incProgress(0.7 / n_files, detail = sprintf("Evaluating [%d/%d]: %s", i, n_files, sim_name))
          
          sim_mat <- read_uploaded_matrix(input$file_multi_sims$datapath[i], input$file_multi_sims$name[i])
          if (i == 1) first_sim_mat <- sim_mat
          
          res_i <- evaluate_simulation_accuracy(
            ref_data = ref_mat,
            sim_data = sim_mat,
            memory_mb = sim_mem,
            elapsed_time = sim_time,
            compute_bivariate = FALSE,
            verbose = FALSE
          )
          
          tbl_i <- res_i$metrics_summary_table
          tbl_i$Method <- sim_name
          if (!"Score" %in% colnames(tbl_i) && "Value" %in% colnames(tbl_i)) {
            tbl_i$Score <- tbl_i$Value
          }
          results_list[[i]] <- tbl_i
        }
        
        incProgress(0.1, detail = "Consolidating all evaluated datasets...")
        combined_df <- do.call(rbind, results_list)
        
        rv$benchmark_df <- combined_df
        rv$methods <- unique(combined_df$Method)
        rv$toy_ref <- ref_mat
        rv$toy_sim <- first_sim_mat
        rv$source_name <- sprintf("Batch Evaluated: %d Simulators against Reference", length(rv$methods))
        
        updateCheckboxGroupInput(session, "sel_bubble_methods", choices = rv$methods, selected = rv$methods)
        showNotification(sprintf("Batch evaluation complete! Evaluated %d simulators simultaneously.", length(rv$methods)), type = "message")
      }, error = function(e) {
        showNotification(paste("Multi-simulator evaluation error:", e$message), type = "error")
      })
    })
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 4 - Load Uploaded Benchmark Object
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_load_uploaded_bench, {
    req(input$file_bench_upload)
    tryCatch({
      ext <- tolower(tools::file_ext(input$file_bench_upload$name))
      if (ext == "rds") {
        obj <- readRDS(input$file_bench_upload$datapath)
        if (is.data.frame(obj)) {
          df <- obj
        } else if (is.list(obj) && !is.null(obj$benchmark_summary_table)) {
          df <- obj$benchmark_summary_table
        } else if (is.list(obj) && !is.null(obj$metrics_summary_table)) {
          df <- obj$metrics_summary_table
        } else if (is.list(obj) && !is.null(obj$summary_table)) {
          df <- obj$summary_table
        } else {
          stop("Unrecognized RDS format. Must contain benchmark summary table.")
        }
      } else if (ext == "csv") {
        df <- utils::read.csv(input$file_bench_upload$datapath, check.names = FALSE)
      } else {
        stop("Unsupported file type. Please upload .rds or .csv.")
      }
      
      if (!"Method" %in% colnames(df) && "Simulator" %in% colnames(df)) df$Method <- df$Simulator
      if (!"Score" %in% colnames(df) && "Value" %in% colnames(df)) df$Score <- df$Value
      
      rv$benchmark_df <- df
      rv$methods <- unique(df$Method)
      rv$source_name <- paste0("Uploaded File: ", input$file_bench_upload$name)
      
      updateCheckboxGroupInput(session, "sel_bubble_methods", choices = rv$methods, selected = rv$methods)
      showNotification("Successfully loaded benchmark results!", type = "message")
    }, error = function(e) {
      showNotification(paste("Upload error:", e$message), type = "error")
    })
  })
  
  # ----------------------------------------------------------------------------
  # Status Banner & Data Preview
  # ----------------------------------------------------------------------------
  output$ui_status_banner <- renderUI({
    if (is.null(rv$benchmark_df)) {
      return(div(class = "alert alert-warning", "No benchmark data loaded yet. Please select an ingestion mode on the left."))
    }
    
    n_methods <- length(unique(rv$benchmark_df$Method))
    n_metrics <- length(unique(rv$benchmark_df$Metric))
    n_records <- nrow(rv$benchmark_df)
    
    div(
      class = "alert alert-success",
      h6(tags$b("Active Benchmark: "), rv$source_name),
      p(sprintf("Total Evaluated Records: %d | Simulators: %d (%s) | Unique Metrics Evaluated: %d",
                n_records, n_methods, paste(rv$methods, collapse = ", "), n_metrics), style = "margin-bottom: 0;")
    )
  })
  
  output$table_active_data_preview <- renderDT({
    req(rv$benchmark_df)
    datatable(
      head(rv$benchmark_df, 50),
      options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE),
      rownames = FALSE,
      class = "compact stripe hover"
    ) %>% formatRound(columns = which(sapply(head(rv$benchmark_df, 50), is.numeric)), digits = 4)
  })
  
  # ----------------------------------------------------------------------------
  # Tab 3: Flagship Comparative Bubble Matrix Controls & Plot
  # ----------------------------------------------------------------------------
  output$ui_bubble_method_picker <- renderUI({
    req(rv$methods)
    checkboxGroupInput(
      "sel_bubble_methods", "Select Simulators to Compare:",
      choices = rv$methods,
      selected = rv$methods
    )
  })
  
  filtered_bubble_data <- reactive({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    
    if (!is.null(input$sel_bubble_methods) && length(input$sel_bubble_methods) > 0) {
      df <- df[df$Method %in% input$sel_bubble_methods, , drop = FALSE]
    }
    
    if (!is.null(input$sel_bubble_cat) && input$sel_bubble_cat != "all") {
      if ("Category" %in% colnames(df)) {
        df <- df[df$Category == input$sel_bubble_cat, , drop = FALSE]
      }
    }
    df
  })
  
  bubble_plot_reactive <- reactive({
    req(filtered_bubble_data())
    df <- filtered_bubble_data()
    req(nrow(df) > 0)
    
    plot_benchmark_bubble_matrix(
      data              = df,
      title             = "scSimEval Studio: Multi-Dimensional Simulation Fidelity Matrix",
      subtitle          = "Direction-aware fidelity scores [0, 1] across curated single-cell evaluation measures",
      base_size         = 9.5,
      bubble_size_range = input$sld_bubble_size,
      show_missing_dots = input$chk_bubble_missing,
      normalize_scores  = input$chk_bubble_norm
    )
  })
  
  output$ui_bubble_plot_container <- renderUI({
    plot_height <- if (!is.null(input$sel_bubble_height)) input$sel_bubble_height else "850px"
    plotOutput("plot_bubble_matrix", height = plot_height)
  })
  
  output$plot_bubble_matrix <- renderPlot({
    bubble_plot_reactive()
  })
  
  # Leaderboard Table Calculation
  leaderboard_reactive <- reactive({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    
    score_col <- if ("Score" %in% colnames(df)) "Score" else if ("Value" %in% colnames(df)) "Value" else NULL
    req(score_col)
    
    w_cat1 <- if (!is.null(input$wt_cat1)) input$wt_cat1 else 1
    w_cat3 <- if (!is.null(input$wt_cat3)) input$wt_cat3 else 1
    w_cat5 <- if (!is.null(input$wt_cat5)) input$wt_cat5 else 1
    w_cat8 <- if (!is.null(input$wt_cat8)) input$wt_cat8 else 1
    
    weights <- rep(1, nrow(df))
    if ("Category" %in% colnames(df)) {
      weights[df$Category == "(I) Distributional Properties"] <- w_cat1
      weights[df$Category == "(III) Cellular Structure & Concordance"] <- w_cat3
      weights[df$Category == "(V) Biological Signal & Downstream Fidelity"] <- w_cat5
      weights[df$Category == "(VIII) Computational Scalability"] <- w_cat8
    }
    
    df$W_Score <- as.numeric(df[[score_col]]) * weights
    
    leaderboard <- aggregate(
      cbind(W_Score, Score = df[[score_col]]) ~ Method,
      data = df,
      FUN = mean,
      na.rm = TRUE
    )
    
    leaderboard$Overall_Rank <- rank(-leaderboard$W_Score, ties.method = "min")
    leaderboard <- leaderboard[order(leaderboard$Overall_Rank), ]
    
    leaderboard$Average_Fidelity <- paste0(round(leaderboard$Score * 100, 1), "%")
    leaderboard$Weighted_Score <- round(leaderboard$W_Score, 3)
    
    leaderboard[, c("Overall_Rank", "Method", "Average_Fidelity", "Weighted_Score")]
  })
  
  output$table_bubble_leaderboard <- renderDT({
    req(leaderboard_reactive())
    datatable(
      leaderboard_reactive(),
      options = list(pageLength = 10, dom = "t"),
      rownames = FALSE,
      class = "compact stripe hover"
    )
  })
  
  # Direct Bubble Plot Downloads
  output$download_bubble_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_bubble_matrix_", Sys.Date(), ".jpeg") },
    content = function(file) {
      export_single_jpeg(file, bubble_plot_reactive(), width = 16, height = 10, dpi = 300)
    }
  )
  output$download_bubble_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_bubble_matrix_", Sys.Date(), ".pdf") },
    content = function(file) {
      grDevices::pdf(file, width = 16, height = 10)
      print(bubble_plot_reactive())
      grDevices::dev.off()
    }
  )
  
  # ----------------------------------------------------------------------------
  # Tab 4: Diagnostic Visualizations (8 Dedicated Sub-Tabs)
  # ----------------------------------------------------------------------------
  
  # 1. Evaluation Summary
  eval_summary_reactive <- reactive({
    req(rv$benchmark_df)
    plot_evaluation_summary(
      data = rv$benchmark_df,
      show_labels = input$chk_sum_labels,
      normalize_scores = input$chk_sum_norm,
      base_size = 11
    )
  })
  output$plot_eval_summary <- renderPlot({ eval_summary_reactive() })
  output$download_sum_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_evaluation_summary_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, eval_summary_reactive(), width = 12, height = 7, dpi = 300) }
  )
  output$download_sum_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_evaluation_summary_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 12, height = 7); print(eval_summary_reactive()); grDevices::dev.off() }
  )
  
  # 2. Distribution QC
  dist_qc_reactive <- reactive({
    req(rv$toy_ref, rv$toy_sim)
    plot_distribution_qc(
      ref_data = rv$toy_ref,
      sim_data = rv$toy_sim,
      layout   = input$sel_dist_layout,
      base_size = 10
    )
  })
  output$plot_dist_qc <- renderPlot({ dist_qc_reactive() })
  output$download_dist_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_distribution_qc_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, dist_qc_reactive(), width = 13, height = 8, dpi = 300) }
  )
  output$download_dist_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_distribution_qc_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 8); print(dist_qc_reactive()); grDevices::dev.off() }
  )
  
  # 3. Scalability Benchmark
  scale_bench_reactive <- reactive({
    req(rv$benchmark_df)
    plot_scalability_benchmark(
      benchmark_data = rv$benchmark_df,
      type = input$sel_scale_type,
      cell_count = input$sld_scale_cells,
      base_size = 11
    )
  })
  output$plot_scale_bench <- renderPlot({ scale_bench_reactive() })
  output$download_scale_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_scalability_benchmark_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, scale_bench_reactive(), width = 13, height = 8, dpi = 300) }
  )
  output$download_scale_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_scalability_benchmark_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 8); print(scale_bench_reactive()); grDevices::dev.off() }
  )
  
  # 4. Metric Boxplots
  metric_box_reactive <- reactive({
    req(rv$benchmark_df)
    plot_metric_boxplots(
      benchmark_data = rv$benchmark_df,
      score_type = input$sel_box_score_type,
      facet_by = input$sel_box_facet,
      base_size = 11
    )
  })
  output$plot_metric_boxes <- renderPlot({ metric_box_reactive() })
  output$download_box_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_boxplots_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, metric_box_reactive(), width = 13, height = 8, dpi = 300) }
  )
  output$download_box_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_boxplots_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 8); print(metric_box_reactive()); grDevices::dev.off() }
  )
  
  # 5. Metric Correlation Heatmap
  metric_heat_reactive <- reactive({
    req(rv$benchmark_df)
    plot_metric_heatmap(
      benchmark_data = rv$benchmark_df,
      cluster_rows = input$chk_heat_cluster_rows,
      cluster_cols = input$chk_heat_cluster_cols,
      base_size = 10
    )
  })
  output$plot_metric_heat <- renderPlot({ metric_heat_reactive() })
  output$download_heat_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_heatmap_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, metric_heat_reactive(), width = 13, height = 8, dpi = 300) }
  )
  output$download_heat_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_heatmap_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 8); print(metric_heat_reactive()); grDevices::dev.off() }
  )
  
  # 6. PCA Ordination
  metric_pca_reactive <- reactive({
    req(rv$benchmark_df)
    plot_metric_pca(
      benchmark_data = rv$benchmark_df,
      panel = input$sel_pca_panel,
      top_n_loadings = input$sld_pca_loadings,
      base_size = 11
    )
  })
  output$plot_metric_pca_out <- renderPlot({ metric_pca_reactive() })
  output$download_pca_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_pca_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, metric_pca_reactive(), width = 13, height = 7, dpi = 300) }
  )
  output$download_pca_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_pca_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 7); print(metric_pca_reactive()); grDevices::dev.off() }
  )
  
  # 7. MDS Metric Space
  metric_mds_reactive <- reactive({
    req(rv$benchmark_df)
    plot_metric_mds(
      benchmark_data = rv$benchmark_df,
      ordination_by = input$sel_mds_by,
      base_size = 11
    )
  })
  output$plot_metric_mds_out <- renderPlot({ metric_mds_reactive() })
  output$download_mds_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_mds_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, metric_mds_reactive(), width = 13, height = 7, dpi = 300) }
  )
  output$download_mds_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_mds_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 7); print(metric_mds_reactive()); grDevices::dev.off() }
  )
  
  # ----------------------------------------------------------------------------
  # Tab 5: Download & Report Center (Excel, CSV, RDS, PDF, and Complete ZIP)
  # ----------------------------------------------------------------------------
  
  # Excel Export (.xlsx)
  output$download_excel <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".xlsx") },
    content = function(file) {
      req(rv$benchmark_df)
      export_excel_workbook(file, rv$benchmark_df, leaderboard_reactive())
    }
  )
  
  # CSV Export (.csv)
  output$download_csv <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".csv") },
    content = function(file) {
      req(rv$benchmark_df)
      utils::write.csv(rv$benchmark_df, file, row.names = FALSE)
    }
  )
  
  # RDS Export (.rds)
  output$download_rds <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".rds") },
    content = function(file) {
      req(rv$benchmark_df)
      saveRDS(list(
        benchmark_summary_table = rv$benchmark_df,
        methods = rv$methods,
        method_rankings = leaderboard_reactive()
      ), file)
    }
  )
  
  # Comprehensive Multi-Page PDF Report
  output$download_all_plots_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_all_plots_report_", Sys.Date(), ".pdf") },
    content = function(file) {
      req(rv$benchmark_df)
      generate_all_plots_pdf(file, rv$benchmark_df, rv$toy_ref, rv$toy_sim)
    }
  )
  
  # Complete Results ZIP Bundle
  output$download_complete_zip <- downloadHandler(
    filename = function() { paste0("scSimEval_complete_benchmark_results_", Sys.Date(), ".zip") },
    content = function(file) {
      req(rv$benchmark_df)
      
      tmp_dir <- file.path(tempdir(), paste0("scSimEval_bundle_", as.integer(Sys.time())))
      dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
      fig_dir <- file.path(tmp_dir, "figures")
      dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
      
      # 1. Excel Workbook (.xlsx)
      export_excel_workbook(file.path(tmp_dir, "scSimEval_benchmark_results.xlsx"), rv$benchmark_df, leaderboard_reactive())
      
      # 2. Master CSV Table (.csv)
      utils::write.csv(rv$benchmark_df, file.path(tmp_dir, "scSimEval_benchmark_results.csv"), row.names = FALSE)
      
      # 3. RDS Object (.rds)
      saveRDS(list(
        benchmark_summary_table = rv$benchmark_df,
        methods = rv$methods,
        method_rankings = leaderboard_reactive()
      ), file.path(tmp_dir, "scSimEval_benchmark_results.rds"))
      
      # 4. Multi-Page PDF Report
      generate_all_plots_pdf(file.path(tmp_dir, "scSimEval_all_plots_report.pdf"), rv$benchmark_df, rv$toy_ref, rv$toy_sim)
      
      # 5. Individual High-Res JPEGs (300 DPI)
      try(export_single_jpeg(file.path(fig_dir, "01_bubble_matrix.jpeg"), bubble_plot_reactive(), width = 16, height = 10, dpi = 300), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "02_evaluation_summary.jpeg"), eval_summary_reactive(), width = 12, height = 7, dpi = 300), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "03_scalability_benchmark.jpeg"), scale_bench_reactive(), width = 12, height = 8, dpi = 300), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "04_metric_boxplots.jpeg"), metric_box_reactive(), width = 12, height = 8, dpi = 300), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "05_metric_heatmap.jpeg"), metric_heat_reactive(), width = 12, height = 8, dpi = 300), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "06_metric_pca.jpeg"), metric_pca_reactive(), width = 12, height = 7, dpi = 300), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "07_metric_mds.jpeg"), metric_mds_reactive(), width = 12, height = 7, dpi = 300), silent = TRUE)
      if (!is.null(rv$toy_ref) && !is.null(rv$toy_sim)) {
        try(export_single_jpeg(file.path(fig_dir, "08_distribution_qc.jpeg"), dist_qc_reactive(), width = 12, height = 8, dpi = 300), silent = TRUE)
      }
      
      # Compress all files into the final .zip
      zip_files <- list.files(tmp_dir, full.names = FALSE, recursive = TRUE)
      zip::zip(file, files = zip_files, root = tmp_dir)
      unlink(tmp_dir, recursive = TRUE)
    }
  )
  
  # Searchable Master Table in Download Center
  output$table_master_export <- renderDT({
    req(rv$benchmark_df)
    datatable(
      rv$benchmark_df,
      filter = "top",
      options = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE,
      class = "compact stripe hover"
    ) %>% formatRound(columns = which(sapply(rv$benchmark_df, is.numeric)), digits = 4)
  })
}

# ==============================================================================
# Shiny App Runner
# ==============================================================================
shinyApp(ui = ui, server = server)
