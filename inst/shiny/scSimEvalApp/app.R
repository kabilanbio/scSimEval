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
      stop("Unsupported RDS format. Please provide a count matrix or data.frame.")
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

# Export multi-sheet Excel workbook
export_excel_workbook <- function(file, benchmark_df, leaderboard_df = NULL) {
  sheets <- list(All_Benchmark_Metrics = benchmark_df)
  if (!is.null(leaderboard_df) && nrow(leaderboard_df) > 0) {
    sheets$Method_Rankings = leaderboard_df
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

# Export High-Res Publication JPEG (600 DPI)
export_single_jpeg <- function(file, plot_obj, width = 14, height = 9, dpi = 600) {
  ggplot2::ggsave(file, plot = plot_obj, device = "jpeg", width = width, height = height, dpi = dpi)
}

# Generate Multi-Page PDF Report
generate_all_plots_pdf <- function(file, benchmark_data, toy_ref = NULL, toy_sim = NULL) {
  grDevices::pdf(file, width = 14, height = 9, onefile = TRUE)
  try(print(plot_benchmark_bubble_matrix(benchmark_data, base_size = 9.5)), silent = TRUE)
  try(print(plot_evaluation_summary(benchmark_data, base_size = 12)), silent = TRUE)
  try(print(plot_scalability_benchmark(benchmark_data, base_size = 12)), silent = TRUE)
  try(print(plot_metric_boxplots(benchmark_data, base_size = 11)), silent = TRUE)
  try(print(plot_metric_heatmap(benchmark_data, base_size = 11)), silent = TRUE)
  try(print(plot_metric_pca(benchmark_data, base_size = 12)), silent = TRUE)
  try(print(plot_metric_mds(benchmark_data, base_size = 12)), silent = TRUE)
  if (!is.null(toy_ref) && !is.null(toy_sim)) {
    try(print(plot_distribution_qc(toy_ref, toy_sim, base_size = 11)), silent = TRUE)
  }
  grDevices::dev.off()
}

# ==============================================================================
# Theme & CSS Styling
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
# UI DEFINITION
# ==============================================================================
ui <- page_navbar(
  title = "scSimEval Studio",
  id = "nav_active",
  theme = app_theme,
  fillable = TRUE,
  
  header = tags$head(
    tags$style(HTML("
      .navbar { box-shadow: 0 2px 8px rgba(0,0,0,0.08); font-weight: 600; }
      .nav-link { font-size: 0.95rem; }
      .stat-card { border-radius: 8px; border-left: 4px solid #1E3A8A; box-shadow: 0 1px 4px rgba(0,0,0,0.05); background: white; padding: 16px; margin-bottom: 15px; }
      .stat-number { font-size: 2.1rem; font-weight: 800; line-height: 1; }
      .stat-label { font-size: 0.8rem; text-transform: uppercase; color: #64748B; font-weight: 600; letter-spacing: 0.5px; margin-top: 5px; }
      .category-pill { display: inline-block; padding: 4px 11px; border-radius: 12px; font-size: 0.8rem; font-weight: 600; color: white; margin: 3px; }
      .hero-box { background: linear-gradient(135deg, #1E3A5F 0%, #243B55 100%); color: white; border-radius: 10px; padding: 24px; margin-bottom: 20px; box-shadow: 0 3px 10px rgba(30,58,95,0.15); }
      .card-header { font-weight: 700; color: #1E293B; background-color: #F8FAFC; border-bottom: 1px solid #E2E8F0; }
      .btn-primary { background-color: #1E3A8A; border-color: #1E3A8A; }
      .btn-primary:hover { background-color: #172554; border-color: #172554; }
      .btn-success { background-color: #0D9488; border-color: #0D9488; }
      .btn-success:hover { background-color: #0F766E; border-color: #0F766E; }
      .guide-step { background: #FFFFFF; border-radius: 8px; border: 1px solid #E2E8F0; padding: 15px; margin-bottom: 12px; }
      .guide-num { display: inline-block; width: 28px; height: 28px; line-height: 28px; border-radius: 50%; background: #1E3A8A; color: white; font-weight: 700; text-align: center; margin-right: 10px; font-size: 0.85rem; }
      
      /* Keep all 7 visual sub-panels in a single non-wrapping row */
      .nav-pills {
        display: flex !important;
        flex-wrap: nowrap !important;
        overflow-x: auto !important;
        white-space: nowrap !important;
        padding-bottom: 6px !important;
        scrollbar-width: thin;
      }
      .nav-pills .nav-item {
        flex: 0 0 auto !important;
      }
      .nav-pills .nav-link {
        font-size: 0.88rem !important;
        padding: 8px 14px !important;
      }
      
      /* Horizontal scroll container for big bubble plot */
      .bubble-scroll-container {
        overflow-x: auto;
        overflow-y: hidden;
        border: 1px solid #E2E8F0;
        border-radius: 8px;
        background: #FFFFFF;
        padding: 12px;
      }
      
      /* Clean scientific inputs */
      .sim-input-card {
        background: #F8FAFC;
        border-left: 4px solid #1B4F72;
        border-radius: 6px;
        padding: 12px;
        margin-bottom: 12px;
      }
    "))
  ),
  
  # ============================================================================
  # TAB 1: HOME
  # ============================================================================
  nav_panel(
    "Home",
    fluidRow(
      column(
        12,
        div(
          class = "hero-box",
          h2("scSimEval: Single-Cell & Multiomics Simulation Benchmarking Studio", style = "font-weight: 800; font-size: 1.85rem; letter-spacing: -0.5px;"),
          p("A unified scientific framework for evaluating and comparing single-cell transcriptomics (scRNA-seq), chromatin accessibility (scATAC-seq), and paired multiomics simulation techniques against empirical biological reference datasets.", style = "font-size: 1.05rem; opacity: 0.95; max-width: 1050px; line-height: 1.5;"),
          hr(style = "border-color: rgba(255,255,255,0.25); margin: 18px 0;"),
          div(
            actionButton("btn_go_data", "1. Data Ingestion & Evaluation", class = "btn btn-outline-light me-2 mb-2", icon = icon("database")),
            actionButton("btn_go_bubble", "2. Comparative Bubble Matrix", class = "btn btn-success me-2 mb-2", icon = icon("chart-pie")),
            actionButton("btn_go_viz", "3. Diagnostic Visualizations", class = "btn btn-info text-white me-2 mb-2", icon = icon("chart-line")),
            actionButton("btn_go_download", "4. Download Results", class = "btn btn-outline-light me-2 mb-2", icon = icon("download")),
            actionButton("btn_go_help", "Documentation & Help", class = "btn btn-outline-light mb-2", icon = icon("book-open"))
          )
        )
      )
    ),
    fluidRow(
      column(3, div(class = "stat-card", style = "border-left-color: #1E3A8A;", div(class = "stat-number", style = "color: #1E3A8A;", "62"), div(class = "stat-label", "Curated Evaluation Measures"))),
      column(3, div(class = "stat-card", style = "border-left-color: #0D9488;", div(class = "stat-number", style = "color: #0D9488;", "8"), div(class = "stat-label", "Canonical Biological Categories"))),
      column(3, div(class = "stat-card", style = "border-left-color: #4F46E5;", div(class = "stat-number", style = "color: #4F46E5;", "3"), div(class = "stat-label", "Supported Data Modalities"))),
      column(3, div(class = "stat-card", style = "border-left-color: #0284C7;", div(class = "stat-number", style = "color: #0284C7;", "0 – 1"), div(class = "stat-label", "Standardized Fidelity Scale")))
    ),
    fluidRow(
      column(
        7,
        card(
          card_header("Eight Evaluation Categories"),
          card_body(
            tags$div(
              style = "margin-bottom: 14px;",
              tags$span(class = "category-pill", style = "background-color: #2563EB;", "(I) Distributional Properties (14 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #0D9488;", "(II) Correlations & Zero-Inflation (6 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #B91C1C;", "(III) Cellular Structure & Concordance (10 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #D97706;", "(IV) Batch Effects & Confounder Mixing (7 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #C2410C;", "(V) Biological Signal & Downstream Fidelity (7 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #7C3AED;", "(VI) Trajectory & Lineage Dynamics (2 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #334155;", "(VII) Cross-Modal Coupling & Modularity (6 metrics)"),
              tags$span(class = "category-pill", style = "background-color: #15803D;", "(VIII) Computational Scalability (2 metrics)")
            ),
            hr(),
            h5("Evaluation Methodology", style = "font-weight: 700; color: #1E293B;"),
            p("The scSimEval framework systematically benchmarks simulation tools across biological, technical, and computational dimensions. By comparing simulated single-cell profiles directly against genuine empirical reference data, it quantifies how faithfully synthetic datasets reproduce true biological properties."),
            p("All raw metrics are transformed through direction-aware standardization to a common [0, 1] scale, allowing seamless multi-metric synthesis, visual matrix comparisons, and objective method rankings.", style = "margin-bottom: 0;")
          )
        )
      ),
      column(
        5,
        card(
          card_header("Standardized Benchmarking Workflow"),
          card_body(
            div(class = "guide-step",
                div(class = "guide-num", "1"),
                tags$b("Data Ingestion: "),
                "Explore pre-computed benchmarks for 6 simulators or upload your own reference empirical counts and simulated datasets (individual or batch)."
            ),
            div(class = "guide-step",
                div(class = "guide-num", "2"),
                tags$b("Comparative Synthesis: "),
                "Examine the 62-metric bubble matrix comparing all simulators side-by-side, along with the automated performance leaderboard."
            ),
            div(class = "guide-step",
                div(class = "guide-num", "3"),
                tags$b("Diagnostic & Publication Deliverables: "),
                "Inspect category summaries, distribution QC, PCA/MDS ordinations, and export high-resolution (600 DPI) figures and Excel workbooks."
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
        width = 400,
        title = "Data Ingestion & Evaluation",
        
        radioButtons(
          "opt_data_mode", "Choose Ingestion Mode:",
          choices = c(
            "Option 1: Explore Demo Benchmark (6 Simulators)" = "demo",
            "Option 2: Evaluate Single Simulator" = "single",
            "Option 3: Evaluate Multiple Simulators (Batch)" = "multi",
            "Option 4: Upload Saved Results (.rds / .csv)" = "upload_bench"
          ),
          selected = "demo"
        ),
        hr(),
        
        # Mode 1: Demo
        conditionalPanel(
          condition = "input.opt_data_mode == 'demo'",
          p("Instantly explore pre-calculated benchmark results for 6 simulation methods (Splatter, scDesign3, SCRIP, SymSim, dyngen, simATAC) across all 62 measures.", style = "font-size: 0.88rem; color: #555;"),
          actionButton("btn_load_demo", "Load Demo Benchmark (6 Simulators)", class = "btn btn-success w-100", icon = icon("play"))
        ),
        
        # Mode 2: Single Simulator
        conditionalPanel(
          condition = "input.opt_data_mode == 'single'",
          h6(tags$b("1. Count Matrices")),
          fileInput("file_single_ref", "Reference Dataset (Real Cells):", accept = c(".rds", ".csv", ".tsv", ".txt")),
          fileInput("file_single_sim", "Simulated Dataset:", accept = c(".rds", ".csv", ".tsv", ".txt")),
          textInput("txt_single_name", "Simulator Method Name:", value = "MySimulator"),
          
          h6(tags$b("2. Computational Scalability (2 Metrics)")),
          fluidRow(
            column(6, numericInput("num_single_time", "Elapsed Time (s):", value = 35.0, min = 0.1, step = 0.5)),
            column(6, numericInput("num_single_mem", "Peak RAM (MB):", value = 820.0, min = 1, step = 10))
          ),
          
          h6(tags$b("3. Optional Annotations")),
          fileInput("file_single_celltypes", "Optional Cell Type Labels (.rds / .csv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
          fileInput("file_single_batch", "Optional Batch Annotations (.rds / .csv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
          
          checkboxInput("chk_append_single", "Append to current benchmark (compare together)", value = TRUE),
          actionButton("btn_run_single_eval", "Run Evaluation", class = "btn btn-primary w-100", icon = icon("calculator"))
        ),
        
        # Mode 3: Multiple Simulators (Batch)
        conditionalPanel(
          condition = "input.opt_data_mode == 'multi'",
          h6(tags$b("1. Reference Biological Dataset")),
          fileInput("file_multi_ref", "Reference Dataset (Real Cells):", accept = c(".rds", ".csv", ".tsv", ".txt")),
          
          h6(tags$b("2. Multiple Simulated Datasets")),
          fileInput("file_multi_sims", "Select Simulated Matrices (2 or more files):", multiple = TRUE, accept = c(".rds", ".csv", ".tsv", ".txt")),
          
          h6(tags$b("3. Optional Annotations (Shared)")),
          fileInput("file_multi_celltypes", "Optional Cell Type Labels (.rds / .csv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
          fileInput("file_multi_batch", "Optional Batch Annotations (.rds / .csv / .txt):", accept = c(".rds", ".csv", ".tsv", ".txt")),
          
          h6(tags$b("4. Configure Simulator Details & Scalability")),
          p("For each selected simulator, specify its name and scalability metrics:", style = "font-size: 0.85rem; color: #555;"),
          uiOutput("ui_multi_sim_scalability_inputs"),
          
          hr(),
          actionButton("btn_run_multi_eval", "Evaluate All Simulators", class = "btn btn-primary w-100", icon = icon("cogs"))
        ),
        
        # Mode 4: Saved Benchmark Upload
        conditionalPanel(
          condition = "input.opt_data_mode == 'upload_bench'",
          p("Upload previously saved evaluation results (.rds or .csv) from scSimEval.", style = "font-size: 0.88rem; color: #555;"),
          fileInput("file_bench_upload", "Select Saved File (.rds or .csv):", accept = c(".rds", ".csv")),
          actionButton("btn_load_uploaded_bench", "Load Saved File", class = "btn btn-info text-white w-100", icon = icon("folder-open"))
        )
      ),
      
      card(
        card_header("Active Benchmark Dataset Status"),
        card_body(
          uiOutput("ui_status_banner"),
          hr(),
          h5("Benchmark Summary Table Preview"),
          p("Displaying evaluated metrics across simulation methods. Column 'Score' represents direction-aware normalized fidelity in [0, 1].", style = "font-size: 0.88rem; color: #7F8C8D;"),
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
        width = 320,
        title = "Display Options",
        
        selectInput(
          "sel_bubble_cat", "Category Filter:",
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
        
        checkboxInput("chk_bubble_missing", "Show Missing Dots (Unmeasured)", value = TRUE),
        
        sliderInput("sld_bubble_width", "Matrix Display Width (px):", min = 800, max = 2600, value = 1400, step = 50),
        p("Tip: Adjust width slider to smoothly scroll horizontally and view all metrics without squishing.", style = "font-size: 0.82rem; color: #666;"),
        
        hr(),
        h6(tags$b("Download This Plot:")),
        downloadButton("download_bubble_jpeg", "Download JPEG (600 DPI)", class = "btn btn-primary w-100 mb-2"),
        downloadButton("download_bubble_pdf", "Download Vector PDF", class = "btn btn-outline-secondary w-100")
      ),
      
      card(
        card_header("Comparative Simulation Fidelity Bubble Matrix (All Datasets)"),
        card_body(
          p("Each column represents a simulator method; each row represents a curated evaluation metric. Bubble size reflects standardized fidelity (larger bubbles = higher fidelity to reference). Color indicates biological category.", style = "font-size: 0.9rem; color: #555;"),
          div(
            class = "bubble-scroll-container",
            uiOutput("ui_bubble_plot_render")
          ),
          hr(),
          h5("Method Ranking Leaderboard (Average Fidelity Score)", style = "font-weight: 700;"),
          DTOutput("table_bubble_leaderboard")
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 4: DIAGNOSTIC VISUALIZATIONS (ALL 7 PANELS IN 1 ROW)
  # ============================================================================
  nav_panel(
    "Visualizations",
    navset_pill(
      id = "viz_subtabs",
      
      # Sub-panel 1: Evaluation Summary
      nav_panel(
        "1. Evaluation Summary",
        card(
          card_header("Overall Category Evaluation Summary (plot_evaluation_summary)"),
          card_body(
            fluidRow(
              column(4, checkboxInput("chk_sum_labels", "Show Score Labels", value = TRUE)),
              column(4, checkboxInput("chk_sum_norm", "Normalize Scores [0, 1]", value = TRUE)),
              column(4,
                     downloadButton("download_sum_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_sum_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_eval_summary", height = "540px")
          )
        )
      ),
      
      # Sub-panel 2: Distribution QC
      nav_panel(
        "2. Distribution QC",
        card(
          card_header("Empirical vs Simulated Distribution Quality (plot_distribution_qc)"),
          card_body(
            fluidRow(
              column(4, selectInput("sel_dist_layout", "QC Layout:", choices = c("Comprehensive" = "comprehensive", "Density Curves Only" = "density"), selected = "comprehensive")),
              column(4, p("Compares expression densities, library sizes, and zero-inflation between real reference and simulated cells.", style = "font-size: 0.85rem; color: #666;")),
              column(4,
                     downloadButton("download_dist_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_dist_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            uiOutput("ui_dist_qc_plot")
          )
        )
      ),
      
      # Sub-panel 3: Scalability Benchmark
      nav_panel(
        "3. Scalability Benchmark",
        card(
          card_header("Computational Scalability Benchmark (plot_scalability_benchmark)"),
          card_body(
            fluidRow(
              column(6,
                     selectInput(
                       "sel_scale_type", "Scalability View:",
                       choices = c(
                         "4-Panel Comprehensive Layout" = "composite",
                         "Runtime Execution Time" = "runtime",
                         "Peak RAM Memory Usage" = "memory",
                         "Runtime vs Memory Trade-Off" = "tradeoff",
                         "Resource Cost Footprint" = "cost"
                       ),
                       selected = "composite"
                     )
              ),
              column(6,
                     downloadButton("download_scale_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_scale_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_scale_bench", height = "620px")
          )
        )
      ),
      
      # Sub-panel 4: Metric Boxplots
      nav_panel(
        "4. Metric Boxplots",
        card(
          card_header("Metric Boxplots & Variance (plot_metric_boxplots)"),
          card_body(
            fluidRow(
              column(12,
                     radioButtons(
                       "opt_box_view_mode", "Select View Mode:",
                       choices = c("View Individual Metric" = "individual", "View by Category Group" = "category"),
                       selected = "individual", inline = TRUE
                     )
              )
            ),
            fluidRow(
              conditionalPanel(
                condition = "input.opt_box_view_mode == 'individual'",
                column(4,
                       selectInput(
                         "sel_box_cat_first", "1. Choose Category First:",
                         choices = c(
                           "(I) Distributional Properties",
                           "(II) Correlations & Zero-Inflation",
                           "(III) Cellular Structure & Concordance",
                           "(IV) Batch Effects & Confounder Mixing",
                           "(V) Biological Signal & Downstream Fidelity",
                           "(VI) Trajectory & Lineage Dynamics",
                           "(VII) Cross-Modal Coupling & Modularity",
                           "(VIII) Computational Scalability"
                         ),
                         selected = "(I) Distributional Properties"
                       )
                ),
                column(4, uiOutput("ui_box_metric_picker")),
                column(4,
                       downloadButton("download_box_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                       downloadButton("download_box_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
                )
              ),
              conditionalPanel(
                condition = "input.opt_box_view_mode == 'category'",
                column(4,
                       selectInput(
                         "sel_box_cat_group", "Choose Category:",
                         choices = c(
                           "All 8 Categories" = "all",
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
                       )
                ),
                column(4, selectInput("sel_box_score_type", "Score Type:", choices = c("Normalized [0, 1]" = "normalized", "Raw Value" = "raw"), selected = "normalized")),
                column(4,
                       downloadButton("download_box_cat_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                       downloadButton("download_box_cat_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
                )
              )
            ),
            hr(),
            plotOutput("plot_metric_boxes", height = "560px")
          )
        )
      ),
      
      # Sub-panel 5: Metric Heatmap
      nav_panel(
        "5. Metric Heatmap",
        card(
          card_header("Metric Correlation & Performance Heatmap (plot_metric_heatmap)"),
          card_body(
            fluidRow(
              column(8, p("Heatmap showing relative performance across all 62 curated metrics and simulation methods.", style = "font-size: 0.9rem; color: #555;")),
              column(4,
                     downloadButton("download_heat_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_heat_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_metric_heat", height = "950px")
          )
        )
      ),
      
      # Sub-panel 6: PCA Ordination
      nav_panel(
        "6. PCA Ordination",
        card(
          card_header("PCA Ordination of Simulators in Performance Space (plot_metric_pca)"),
          card_body(
            fluidRow(
              column(4,
                     selectInput(
                       "sel_pca_cat", "Category Filter (6 Categories):",
                       choices = c(
                         "All Categories Combined" = "all",
                         "(I) Distributional Properties" = "(I) Distributional Properties",
                         "(II) Correlations & Zero-Inflation" = "(II) Correlations & Zero-Inflation",
                         "(III) Cellular Structure & Concordance" = "(III) Cellular Structure & Concordance",
                         "(IV) Batch Effects & Confounder Mixing" = "(IV) Batch Effects & Confounder Mixing",
                         "(V) Biological Signal & Downstream Fidelity" = "(V) Biological Signal & Downstream Fidelity",
                         "(VII) Cross-Modal Coupling & Modularity" = "(VII) Cross-Modal Coupling & Modularity"
                       ),
                       selected = "all"
                     )
              ),
              column(4, selectInput("sel_pca_panel", "Panel View:", choices = c("Both (Biplot + Loadings)" = "both", "Simulators Only" = "methods", "Loadings Only" = "loadings"), selected = "both")),
              column(4,
                     downloadButton("download_pca_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_pca_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_metric_pca_out", height = "560px")
          )
        )
      ),
      
      # Sub-panel 7: MDS Metric Space
      nav_panel(
        "7. MDS Metric Space",
        card(
          card_header("Multi-Dimensional Scaling (MDS) Ordination (plot_metric_mds)"),
          card_body(
            fluidRow(
              column(4,
                     selectInput(
                       "sel_mds_cat", "Category Filter (6 Categories):",
                       choices = c(
                         "All Categories Combined" = "all",
                         "(I) Distributional Properties" = "(I) Distributional Properties",
                         "(II) Correlations & Zero-Inflation" = "(II) Correlations & Zero-Inflation",
                         "(III) Cellular Structure & Concordance" = "(III) Cellular Structure & Concordance",
                         "(IV) Batch Effects & Confounder Mixing" = "(IV) Batch Effects & Confounder Mixing",
                         "(V) Biological Signal & Downstream Fidelity" = "(V) Biological Signal & Downstream Fidelity",
                         "(VII) Cross-Modal Coupling & Modularity" = "(VII) Cross-Modal Coupling & Modularity"
                       ),
                       selected = "all"
                     )
              ),
              column(4, selectInput("sel_mds_by", "MDS Target:", choices = c("By Simulators" = "simulators", "By Metric Summaries" = "summaries"), selected = "simulators")),
              column(4,
                     downloadButton("download_mds_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2"),
                     downloadButton("download_mds_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary")
              )
            ),
            hr(),
            plotOutput("plot_metric_mds_out", height = "560px")
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 5: DOWNLOAD RESULTS (BEGINNER-FRIENDLY)
  # ============================================================================
  nav_panel(
    "Download Results",
    fluidRow(
      column(
        4,
        card(
          card_header("1. All-in-One Benchmark Archive (.zip)"),
          card_body(
            p("Download everything at once in a single convenient zip file:", style = "font-size: 0.92rem;"),
            tags$ul(
              tags$li(tags$b("Excel Workbook (.xlsx): "), "Full results with method rankings."),
              tags$li(tags$b("Master Table (.csv): "), "All 62 metrics in tidy format."),
              tags$li(tags$b("R Object (.rds): "), "For downstream R analysis."),
              tags$li(tags$b("Multi-Page PDF Report: "), "All 7 figures compiled."),
              tags$li(tags$b("Figure Images: "), "Individual 600 DPI publication JPEGs.")
            ),
            hr(),
            downloadButton("download_complete_zip", "Download Complete Results (.zip)", class = "btn btn-success w-100 py-2", icon = icon("file-zipper"))
          )
        )
      ),
      column(
        4,
        card(
          card_header("2. Spreadsheets & Data Files"),
          card_body(
            p("Open and analyze your evaluation scores in Microsoft Excel, Google Sheets, or R:", style = "font-size: 0.92rem;"),
            downloadButton("download_excel", "Download Excel File (.xlsx)", class = "btn btn-primary w-100 mb-2", icon = icon("file-excel")),
            downloadButton("download_csv", "Download CSV Table (.csv)", class = "btn btn-outline-primary w-100 mb-2", icon = icon("file-csv")),
            downloadButton("download_rds", "Download R Data File (.rds)", class = "btn btn-outline-secondary w-100", icon = icon("code"))
          )
        )
      ),
      column(
        4,
        card(
          card_header("3. Complete Multi-Page PDF Report"),
          card_body(
            p("Download all evaluation figures compiled into a single high-quality PDF report:", style = "font-size: 0.92rem;"),
            tags$ol(
              tags$li("Comparative Bubble Matrix"),
              tags$li("Evaluation Summary"),
              tags$li("Scalability Benchmark"),
              tags$li("Metric Boxplots"),
              tags$li("Metric Heatmap"),
              tags$li("PCA Ordination"),
              tags$li("MDS Metric Space"),
              tags$li("Distribution QC Curves")
            ),
            hr(),
            downloadButton("download_all_plots_pdf", "Download All Figures (.pdf)", class = "btn btn-info text-white w-100 py-2", icon = icon("file-pdf"))
          )
        )
      )
    ),
    fluidRow(
      column(
        12,
        card(
          card_header("Interactive Benchmark Data Table"),
          card_body(
            p("Type in the search boxes below to automatically filter by Method Name, Category, or Metric:", style = "font-size: 0.9rem; color: #555;"),
            DTOutput("table_master_export")
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 6: HELP & GETTING STARTED
  # ============================================================================
  nav_panel(
    "Help & Getting Started",
    fluidRow(
      column(
        12,
        card(
          card_header("Getting Started with scSimEval"),
          card_body(
            h4("1. Overview", style = "font-weight: 700; color: #1B4F72;"),
            p("Computer simulations of single-cell technologies (scRNA-seq, scATAC-seq, and paired multiomics) are widely used to test bioinformatics pipelines, benchmark statistical tools, and evaluate experimental designs. A central question is always: ",
              tags$i("how realistic is the simulated data compared to genuine biological experiments?")),
            p(tags$b("scSimEval"), " provides ", tags$b("62 evaluation measures organized into 8 easy-to-understand categories"), ". It evaluates simulation techniques directly against real empirical datasets without requiring artificial ground-truth labels."),
            hr(),
            
            h4("2. The Eight Evaluation Categories", style = "font-weight: 700; color: #1B4F72;"),
            tags$table(
              class = "table table-bordered table-striped",
              tags$thead(
                tags$tr(
                  tags$th("Category"),
                  tags$th("Metrics Count"),
                  tags$th("What It Evaluates"),
                  tags$th("Key Measures Included")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td(tags$b("(I) Distributional Properties")),
                  tags$td("14 metrics"),
                  tags$td("Statistical distance between real and simulated expression distributions at both cell and feature levels."),
                  tags$td("KS distance, Wasserstein distance, MAD, MAE, RMSE, Bhattacharyya distance, Jaccard distance")
                ),
                tags$tr(
                  tags$td(tags$b("(II) Correlations & Zero-Inflation")),
                  tags$td("6 metrics"),
                  tags$td("Gene-gene co-expression, cell-cell correlations, and dropout patterns."),
                  tags$td("Gene correlation difference, cell correlation difference, zero fraction divergence")
                ),
                tags$tr(
                  tags$td(tags$b("(III) Cellular Structure & Concordance")),
                  tags$td("10 metrics"),
                  tags$td("How faithfully cell types, clustering boundaries, and manifold geometry are preserved."),
                  tags$td("Silhouette width, Adjusted Rand Index (ARI), Normalized Mutual Information (NMI), Neighborhood Purity")
                ),
                tags$tr(
                  tags$td(tags$b("(IV) Batch Effects & Confounder Mixing")),
                  tags$td("7 metrics"),
                  tags$td("Evaluation of technical batch variation and biological cell-type mixing."),
                  tags$td("kBET rejection rate, Batch LISI, Cell-type LISI, CMS score, Mixing metric")
                ),
                tags$tr(
                  tags$td(tags$b("(V) Biological Signal & Downstream Fidelity")),
                  tags$td("7 metrics"),
                  tags$td("Preservation of biological marker genes and differential expression (DEG) rankings."),
                  tags$td("DEG overlap, Jaccard index, Spearman rank correlation of logFC, F1-score")
                ),
                tags$tr(
                  tags$td(tags$b("(VI) Trajectory & Lineage Dynamics")),
                  tags$td("2 metrics"),
                  tags$td("Preservation of continuous developmental pathways and pseudotime progression."),
                  tags$td("Trajectory topology concordance, Pseudotime Spearman correlation")
                ),
                tags$tr(
                  tags$td(tags$b("(VII) Cross-Modal Coupling & Modularity")),
                  tags$td("6 metrics"),
                  tags$td("Coordination between paired modalities (e.g. gene expression and chromatin accessibility)."),
                  tags$td("Cross-modal correlation, Modality concordance, Paired cell distance")
                ),
                tags$tr(
                  tags$td(tags$b("(VIII) Computational Scalability")),
                  tags$td("2 metrics"),
                  tags$td("Computational efficiency and resource usage."),
                  tags$td("Elapsed runtime (seconds), Peak memory usage (MB)")
                )
              )
            ),
            hr(),
            
            h4("3. Two-Step Score Normalization Pipeline", style = "font-weight: 700; color: #1B4F72;"),
            p("In single-cell benchmarking, different metrics have different units and directions. For example, runtime is in seconds, peak memory is in megabytes, statistical distances are near zero, and clustering accuracy ranges between -1 and 1. For some metrics, smaller values are better (error, runtime), while for others, larger values are better (correlation, ARI)."),
            p("To make fair comparisons, ", tags$b("scSimEval"), " applies a standardized two-step normalization:"),
            tags$ol(
              tags$li(tags$b("Direction Inversion: "), "Metrics where lower values indicate better results are inverted so higher scores always indicate superior performance: Inverted = Max - Value."),
              tags$li(tags$b("Min-Max Scaling [0.00, 1.00]: "), "Scores are scaled between 0 (worst performer) and 1 (best performer): Score = (Value - Min) / (Max - Min).")
            ),
            hr(),
            
            h4("4. How to Run in R (Code Examples)", style = "font-weight: 700; color: #1B4F72;"),
            p("You can execute the exact same benchmarking workflows directly in R:"),
            tags$pre(
              tags$code(
                "# 1. Unimodal scRNA-seq Simulation Accuracy:\n",
                "library(scSimEval)\n",
                "results <- evaluate_simulation_accuracy(\n",
                "  ref_data     = real_counts_matrix,\n",
                "  sim_data     = simulated_counts_matrix,\n",
                "  elapsed_time = 45.2,   # seconds\n",
                "  memory_mb    = 850     # peak RAM in MB\n",
                ")\n\n",
                "# 2. Launch this Interactive Shiny Studio:\n",
                "launch_scSimEval_app()\n\n",
                "# 3. Plot the Comparative Bubble Matrix:\n",
                "plot_benchmark_bubble_matrix(results$metrics_summary_table)\n"
              )
            ),
            hr(),
            
            h4("5. Frequently Asked Questions (FAQ)", style = "font-weight: 700; color: #1B4F72;"),
            tags$ul(
              tags$li(tags$b("What file formats are supported? "), "You can upload .rds (matrices or data frames), .csv, .tsv, or .txt files."),
              tags$li(tags$b("Do I need artificial ground truth labels? "), "No. scSimEval evaluates how well simulated data reproduce genuine biological reference datasets across statistical, cellular, and molecular properties."),
              tags$li(tags$b("How should I measure runtime and memory? "), "Record the wall-clock execution time (seconds) and the peak resident memory (MB) consumed by your simulator, then enter them in the Scalability inputs.")
            )
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
    source_name = if (!is.null(initial_demo)) "Demo Benchmark (Splatter, scDesign3, SCRIP, SymSim, dyngen, simATAC)" else "No Data Loaded"
  )
  
  # Navigation triggers
  observeEvent(input$btn_go_data, { nav_select("nav_active", "Data Hub") })
  observeEvent(input$btn_go_bubble, { nav_select("nav_active", "Comparative Bubble Matrix") })
  observeEvent(input$btn_go_viz, { nav_select("nav_active", "Visualizations") })
  observeEvent(input$btn_go_download, { nav_select("nav_active", "Download Results") })
  observeEvent(input$btn_go_help, { nav_select("nav_active", "Help & Getting Started") })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 1 - Load Demo Benchmark
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_load_demo, {
    if (!is.null(initial_demo)) {
      rv$benchmark_df <- initial_demo$benchmark_summary_table
      rv$methods <- initial_demo$methods
      rv$toy_ref <- initial_demo$toy_data$ref
      rv$toy_sim <- initial_demo$toy_data$sim
      rv$source_name <- "Demo Benchmark (Splatter, scDesign3, SCRIP, SymSim, dyngen, simATAC)"
      
      updateCheckboxGroupInput(session, "sel_bubble_methods", choices = rv$methods, selected = rv$methods)
      showNotification("Demo benchmark loaded successfully!", type = "message")
    } else {
      showNotification("Demo benchmark file not found on disk.", type = "warning")
    }
  })
  
  # ----------------------------------------------------------------------------
  # Data Hub: Dynamic Inputs for Multiple Simulators Upload (Mode 3)
  # ----------------------------------------------------------------------------
  output$ui_multi_sim_scalability_inputs <- renderUI({
    req(input$file_multi_sims)
    n_files <- nrow(input$file_multi_sims)
    if (n_files == 0) return(NULL)
    
    inputs_list <- lapply(seq_len(n_files), function(i) {
      fname <- input$file_multi_sims$name[i]
      default_name <- tools::file_path_sans_ext(fname)
      
      div(
        class = "sim-input-card",
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
    
    withProgress(message = "Evaluating simulator...", detail = "Loading count matrices", value = 0.2, {
      tryCatch({
        ref_mat <- read_uploaded_matrix(input$file_single_ref$datapath, input$file_single_ref$name)
        sim_mat <- read_uploaded_matrix(input$file_single_sim$datapath, input$file_single_sim$name)
        
        incProgress(0.3, detail = "Computing 62 evaluation metrics...")
        
        sim_name <- trimws(input$txt_single_name)
        if (sim_name == "") sim_name <- "Simulator"
        
        elapsed_sec <- as.numeric(input$num_single_time)
        peak_ram <- as.numeric(input$num_single_mem)
        
        cell_types_vec <- read_uploaded_labels(input$file_single_celltypes$datapath, input$file_single_celltypes$name)
        batch_vec <- read_uploaded_labels(input$file_single_batch$datapath, input$file_single_batch$name)
        
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
        
        if (!"Score" %in% colnames(tbl) && "Value" %in% colnames(tbl)) {
          tbl$Score <- tbl$Value
        }
        
        if (isTRUE(input$chk_append_single) && !is.null(rv$benchmark_df)) {
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
    
    withProgress(message = "Batch Evaluation", value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Loading biological reference matrix...")
        ref_mat <- read_uploaded_matrix(input$file_multi_ref$datapath, input$file_multi_ref$name)
        
        cell_types_vec <- read_uploaded_labels(input$file_multi_celltypes$datapath, input$file_multi_celltypes$name)
        batch_vec <- read_uploaded_labels(input$file_multi_batch$datapath, input$file_multi_batch$name)
        
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
  # Data Hub: Mode 4 - Load Saved Benchmark File
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
      showNotification("Benchmark results loaded successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Upload error:", e$message), type = "error")
    })
  })
  
  # ----------------------------------------------------------------------------
  # Status Banner & Data Preview
  # ----------------------------------------------------------------------------
  output$ui_status_banner <- renderUI({
    if (is.null(rv$benchmark_df)) {
      return(div(class = "alert alert-warning", "No benchmark data loaded yet. Please select an option on the left."))
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
  # Tab 3: Comparative Bubble Matrix
  # ----------------------------------------------------------------------------
  output$ui_bubble_method_picker <- renderUI({
    req(rv$methods)
    checkboxGroupInput(
      "sel_bubble_methods", "Select Simulators:",
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
      subtitle          = "Standardized Direction-Aware Fidelity Scores [0, 1] Across 8 Canonical Evaluation Categories",
      base_size         = 10.5,
      show_missing_dots = input$chk_bubble_missing,
      normalize_scores  = TRUE
    )
  })
  
  output$ui_bubble_plot_render <- renderUI({
    w <- if (!is.null(input$sld_bubble_width)) paste0(input$sld_bubble_width, "px") else "1400px"
    plotOutput("plot_bubble_matrix", width = w, height = "850px")
  })
  
  output$plot_bubble_matrix <- renderPlot({
    bubble_plot_reactive()
  })
  
  # Method Ranking Leaderboard
  leaderboard_reactive <- reactive({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    score_col <- if ("Score" %in% colnames(df)) "Score" else if ("Value" %in% colnames(df)) "Value" else NULL
    req(score_col)
    
    leaderboard <- aggregate(
      df[[score_col]],
      by = list(Method = df$Method),
      FUN = mean,
      na.rm = TRUE
    )
    colnames(leaderboard)[2] <- "Score"
    
    leaderboard$Overall_Rank <- rank(-leaderboard$Score, ties.method = "min")
    leaderboard <- leaderboard[order(leaderboard$Overall_Rank), ]
    leaderboard$Average_Fidelity <- paste0(round(leaderboard$Score * 100, 1), "%")
    leaderboard$Fidelity_Score <- round(leaderboard$Score, 4)
    
    leaderboard[, c("Overall_Rank", "Method", "Average_Fidelity", "Fidelity_Score")]
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
  
  output$download_bubble_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_bubble_matrix_", Sys.Date(), ".jpeg") },
    content = function(file) {
      export_single_jpeg(file, bubble_plot_reactive(), width = 16, height = 10, dpi = 600)
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
  # Tab 4: Diagnostic Visualizations (7 Panels in 1 Row)
  # ----------------------------------------------------------------------------
  
  # 1. Evaluation Summary
  eval_summary_reactive <- reactive({
    req(rv$benchmark_df)
    plot_evaluation_summary(
      data = rv$benchmark_df,
      show_labels = input$chk_sum_labels,
      normalize_scores = input$chk_sum_norm,
      base_size = 14
    )
  })
  output$plot_eval_summary <- renderPlot({ eval_summary_reactive() })
  output$download_sum_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_evaluation_summary_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, eval_summary_reactive(), width = 13, height = 7.5, dpi = 600) }
  )
  output$download_sum_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_evaluation_summary_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 7.5); print(eval_summary_reactive()); grDevices::dev.off() }
  )
  
  # 2. Distribution QC
  output$ui_dist_qc_plot <- renderUI({
    plot_h <- if (identical(input$sel_dist_layout, "comprehensive")) "850px" else "550px"
    plotOutput("plot_dist_qc", height = plot_h)
  })
  dist_qc_reactive <- reactive({
    req(rv$toy_ref, rv$toy_sim)
    plot_distribution_qc(
      ref_data = rv$toy_ref,
      sim_data = rv$toy_sim,
      layout   = input$sel_dist_layout,
      base_size = 13
    )
  })
  output$plot_dist_qc <- renderPlot({ dist_qc_reactive() })
  output$download_dist_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_distribution_qc_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, dist_qc_reactive(), width = 14, height = 9, dpi = 600) }
  )
  output$download_dist_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_distribution_qc_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 14, height = 9); print(dist_qc_reactive()); grDevices::dev.off() }
  )
  
  # 3. Scalability Benchmark
  scale_bench_reactive <- reactive({
    req(rv$benchmark_df)
    plot_scalability_benchmark(
      benchmark_data = rv$benchmark_df,
      type = input$sel_scale_type,
      base_size = 13
    )
  })
  output$plot_scale_bench <- renderPlot({ scale_bench_reactive() })
  output$download_scale_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_scalability_benchmark_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, scale_bench_reactive(), width = 13, height = 8, dpi = 600) }
  )
  output$download_scale_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_scalability_benchmark_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 8); print(scale_bench_reactive()); grDevices::dev.off() }
  )
  
  # 4. Metric Boxplots
  output$ui_box_metric_picker <- renderUI({
    req(rv$benchmark_df, input$sel_box_cat_first)
    sub_df <- rv$benchmark_df[rv$benchmark_df$Category == input$sel_box_cat_first, , drop = FALSE]
    avail_metrics <- sort(unique(sub_df$Metric))
    selectInput("sel_box_metric_single", "2. Choose Metric Name:", choices = avail_metrics, selected = avail_metrics[1])
  })
  
  metric_box_reactive <- reactive({
    req(rv$benchmark_df)
    if (identical(input$opt_box_view_mode, "individual")) {
      req(input$sel_box_metric_single)
      plot_metric_boxplots(
        benchmark_data = rv$benchmark_df,
        metrics = input$sel_box_metric_single,
        score_type = "normalized",
        base_size = 12
      )
    } else {
      cat_filter <- if (identical(input$sel_box_cat_group, "all")) NULL else input$sel_box_cat_group
      plot_metric_boxplots(
        benchmark_data = rv$benchmark_df,
        categories = cat_filter,
        score_type = input$sel_box_score_type,
        facet_by = if (is.null(cat_filter)) "category" else "metric",
        base_size = 11
      )
    }
  })
  output$plot_metric_boxes <- renderPlot({ metric_box_reactive() })
  output$download_box_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_boxplot_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, metric_box_reactive(), width = 13, height = 7.5, dpi = 600) }
  )
  output$download_box_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_boxplot_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 7.5); print(metric_box_reactive()); grDevices::dev.off() }
  )
  output$download_box_cat_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_category_boxplots_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, metric_box_reactive(), width = 13, height = 8, dpi = 600) }
  )
  output$download_box_cat_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_category_boxplots_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 8); print(metric_box_reactive()); grDevices::dev.off() }
  )
  
  # 5. Metric Correlation Heatmap (Increased height, no clustering options)
  metric_heat_reactive <- reactive({
    req(rv$benchmark_df)
    plot_metric_heatmap(
      benchmark_data = rv$benchmark_df,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      base_size = 12
    )
  })
  output$plot_metric_heat <- renderPlot({ metric_heat_reactive() })
  output$download_heat_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_heatmap_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, metric_heat_reactive(), width = 14, height = 12, dpi = 600) }
  )
  output$download_heat_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_heatmap_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 14, height = 12); print(metric_heat_reactive()); grDevices::dev.off() }
  )
  
  # 6. PCA Ordination (6 Category-wise options)
  metric_pca_reactive <- reactive({
    req(rv$benchmark_df)
    cat_sel <- if (identical(input$sel_pca_cat, "all")) NULL else input$sel_pca_cat
    plot_metric_pca(
      benchmark_data = rv$benchmark_df,
      category = cat_sel,
      panel = input$sel_pca_panel,
      base_size = 13
    )
  })
  output$plot_metric_pca_out <- renderPlot({ metric_pca_reactive() })
  output$download_pca_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_pca_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, metric_pca_reactive(), width = 13, height = 7.5, dpi = 600) }
  )
  output$download_pca_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_pca_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 7.5); print(metric_pca_reactive()); grDevices::dev.off() }
  )
  
  # 7. MDS Metric Space (6 Category-wise options)
  metric_mds_reactive <- reactive({
    req(rv$benchmark_df)
    cat_sel <- if (identical(input$sel_mds_cat, "all")) NULL else input$sel_mds_cat
    plot_metric_mds(
      benchmark_data = rv$benchmark_df,
      category = cat_sel,
      ordination_by = input$sel_mds_by,
      base_size = 13
    )
  })
  output$plot_metric_mds_out <- renderPlot({ metric_mds_reactive() })
  output$download_mds_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_mds_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, metric_mds_reactive(), width = 13, height = 7.5, dpi = 600) }
  )
  output$download_mds_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_mds_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 7.5); print(metric_mds_reactive()); grDevices::dev.off() }
  )
  
  # ----------------------------------------------------------------------------
  # Tab 5: Download Results (Excel, CSV, RDS, PDF, and Complete ZIP)
  # ----------------------------------------------------------------------------
  
  output$download_excel <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".xlsx") },
    content = function(file) {
      req(rv$benchmark_df)
      export_excel_workbook(file, rv$benchmark_df, leaderboard_reactive())
    }
  )
  
  output$download_csv <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".csv") },
    content = function(file) {
      req(rv$benchmark_df)
      utils::write.csv(rv$benchmark_df, file, row.names = FALSE)
    }
  )
  
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
  
  output$download_all_plots_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_all_plots_report_", Sys.Date(), ".pdf") },
    content = function(file) {
      req(rv$benchmark_df)
      generate_all_plots_pdf(file, rv$benchmark_df, rv$toy_ref, rv$toy_sim)
    }
  )
  
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
      
      # 5. Publication-Ready JPEGs at 600 DPI
      try(export_single_jpeg(file.path(fig_dir, "01_bubble_matrix.jpeg"), bubble_plot_reactive(), width = 16, height = 10, dpi = 600), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "02_evaluation_summary.jpeg"), eval_summary_reactive(), width = 13, height = 7.5, dpi = 600), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "03_scalability_benchmark.jpeg"), scale_bench_reactive(), width = 13, height = 8, dpi = 600), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "04_metric_boxplots.jpeg"), metric_box_reactive(), width = 13, height = 7.5, dpi = 600), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "05_metric_heatmap.jpeg"), metric_heat_reactive(), width = 14, height = 12, dpi = 600), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "06_metric_pca.jpeg"), metric_pca_reactive(), width = 13, height = 7.5, dpi = 600), silent = TRUE)
      try(export_single_jpeg(file.path(fig_dir, "07_metric_mds.jpeg"), metric_mds_reactive(), width = 13, height = 7.5, dpi = 600), silent = TRUE)
      if (!is.null(rv$toy_ref) && !is.null(rv$toy_sim)) {
        try(export_single_jpeg(file.path(fig_dir, "08_distribution_qc.jpeg"), dist_qc_reactive(), width = 14, height = 9, dpi = 600), silent = TRUE)
      }
      
      zip_files <- list.files(tmp_dir, full.names = FALSE, recursive = TRUE)
      zip::zip(file, files = zip_files, root = tmp_dir)
      unlink(tmp_dir, recursive = TRUE)
    }
  )
  
  # Searchable Master Table with Factor Filter Dropdowns
  output$table_master_export <- renderDT({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    
    # Convert text columns to factors for automatic dropdown filtering
    if ("Method" %in% colnames(df)) df$Method <- as.factor(df$Method)
    if ("Category" %in% colnames(df)) df$Category <- as.factor(df$Category)
    if ("Metric" %in% colnames(df)) df$Metric <- as.factor(df$Metric)
    
    datatable(
      df,
      filter = list(position = "top", clear = FALSE),
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        autoWidth = TRUE,
        searchHighlight = TRUE
      ),
      rownames = FALSE,
      class = "compact stripe hover"
    ) %>% formatRound(columns = which(sapply(df, is.numeric)), digits = 4)
  })
}

# ==============================================================================
# Run Shiny App
# ==============================================================================
shinyApp(ui = ui, server = server)
