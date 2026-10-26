# inst/shiny/scSimEvalApp/app.R
# Simulation Evaluation Studio for scRNA-seq, scATAC-seq, and Single-Cell Multiomics
# Powered by scSimEval (62 Curated Benchmark Measures)

library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(Matrix)
library(scSimEval)

# If running in local checkout, source updated visualizations to guarantee latest bugfixes
for (p in c("10_visualizations.R", "R/10_visualizations.R", "../../R/10_visualizations.R", "../../../R/10_visualizations.R")) {
  if (file.exists(p)) {
    try(source(p, local = FALSE), silent = TRUE)
    break
  }
}

# Set max upload size to 10 GB for large single-cell datasets
options(shiny.maxRequestSize = 10 * 1024^3)

# Detect available CPU hardware cores for parallel benchmarking
max_sys_cores <- tryCatch({
  mc <- parallel::detectCores(logical = TRUE)
  if (is.na(mc) || mc < 1) 1L else as.integer(mc)
}, error = function(e) 1L)
default_n_cores <- max(1L, min(4L, max_sys_cores - 1L))


# Register static resource directory for documentation figures
fig_dirs <- c(
  system.file("shiny", "scSimEvalApp", "www", package = "scSimEval"),
  system.file("www", package = "scSimEval"),
  file.path("inst", "shiny", "scSimEvalApp", "www"),
  file.path("www"),
  file.path("man", "figures"),
  file.path("..", "man", "figures"),
  file.path("..", "..", "man", "figures")
)
for (fd in fig_dirs) {
  if (nzchar(fd) && dir.exists(fd)) {
    try(shiny::addResourcePath("scfigures", normalizePath(fd, winslash = "/", mustWork = FALSE)), silent = TRUE)
    break
  }
}

# Load demo benchmark data if present
demo_data_path <- system.file("shiny", "scSimEvalApp", "data", "demo_benchmark_data.rds", package = "scSimEval")
if (demo_data_path == "" || !file.exists(demo_data_path)) {
  demo_data_path <- file.path("data", "demo_benchmark_data.rds")
}
initial_demo <- if (file.exists(demo_data_path)) readRDS(demo_data_path) else NULL

# Source benchmark suite evaluator if running locally
for (p in c("09_evaluate_benchmark.R", "R/09_evaluate_benchmark.R", "../../R/09_evaluate_benchmark.R", "../../../R/09_evaluate_benchmark.R")) {
  if (file.exists(p)) {
    try(source(p, local = FALSE), silent = TRUE)
    break
  }
}

# Standardize Category names to canonical 8 package categories
standardize_benchmark_categories <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(df)
  legacy_map <- c(
    "Distributional Properties"               = "(I) Distributional Properties",
    "Distribution"                            = "(I) Distributional Properties",
    "Correlation & Dependencies"              = "(II) Correlations & Zero-Inflation",
    "Correlations & Zero-Inflation"           = "(II) Correlations & Zero-Inflation",
    "Correlation"                             = "(II) Correlations & Zero-Inflation",
    "Cellular Structure & Mixing"             = "(III) Cellular Structure & Concordance",
    "Cellular Structure & Concordance"        = "(III) Cellular Structure & Concordance",
    "Cell Structure"                          = "(III) Cellular Structure & Concordance",
    "Batch Effects & Confounder Mixing"       = "(IV) Batch Effects & Confounder Mixing",
    "Batch Mixing"                            = "(IV) Batch Effects & Confounder Mixing",
    "Biological Signal & Downstream"          = "(V) Biological Signal & Downstream Fidelity",
    "Biological Signal & Downstream Fidelity" = "(V) Biological Signal & Downstream Fidelity",
    "Bio-Signal & DE"                         = "(V) Biological Signal & Downstream Fidelity",
    "Trajectory Dynamics"                     = "(VI) Trajectory & Lineage Dynamics",
    "Trajectory & Lineage Dynamics"           = "(VI) Trajectory & Lineage Dynamics",
    "Trajectory"                              = "(VI) Trajectory & Lineage Dynamics",
    "Traj."                                   = "(VI) Trajectory & Lineage Dynamics",
    "Cross-Modal Relationships"               = "(VII) Cross-Modal Coupling & Modularity",
    "Cross-Modal Coupling & Modularity"       = "(VII) Cross-Modal Coupling & Modularity",
    "Cross-Modal"                             = "(VII) Cross-Modal Coupling & Modularity",
    "Computational Scalability"               = "(VIII) Computational Scalability",
    "Scalability"                             = "(VIII) Computational Scalability"
  )
  if ("Category" %in% colnames(df)) {
    df$Category <- ifelse(df$Category %in% names(legacy_map),
                          legacy_map[df$Category], df$Category)
  }
  if ("Metric" %in% colnames(df)) {
    mc <- tryCatch(scSimEval:::.METRIC_CATEGORY_MAP[df$Metric], error = function(e) rep(NA_character_, nrow(df)))
    batch_metrics <- c("batch_silhouette", "shannon_entropy", "pcr_r2", "cms", "isi", "seurat_mixing_metric", "ldf_diff")
    zero_metrics  <- c("midpoint_error", "slope_error", "r2_discrepancy", "excess_zero_weight_discrepancy", "dropout_curve_midpoint_x0", "dropout_curve_slope_k", "mean_variance_r2", "bcv_discrepancy")
    
    needs_reassign <- !is.na(mc) & (
      is.na(df$Category) | df$Category == "" | df$Category == "Uncategorized" |
      (!is.na(df$Category) & df$Metric %in% batch_metrics & df$Category != "(IV) Batch Effects & Confounder Mixing") |
      (!is.na(df$Category) & df$Metric %in% zero_metrics & df$Category != "(II) Correlations & Zero-Inflation")
    )
    needs_reassign[is.na(needs_reassign)] <- FALSE
    if (any(needs_reassign)) {
      df$Category[needs_reassign] <- mc[needs_reassign]
    }
  }
  df
}

.HIGHER_IS_BETTER_METRICS_LIST <- c(
  "OV", "mean_variance_r2", "cell_pearson_cor", "feature_pearson_cor", 
  "silhouette_sim", "dunn_sim", "calinski_harabasz_sim", "ari", "nmi", "ami", "v_measure", "neighborhood_purity", "Generative_Precision", 
  "cms", "isi", "seurat_mixing_metric", "shannon_entropy", 
  "SimBench_DE_Fidelity_Score", "Log2FC_Pearson_Corr", "Log2FC_Spearman_Corr", "Top_DEG_Jaccard_Overlap", 
  "Distribution_Score", "Classifier_Accuracy", "Classifier_Macro_F1", "Classifier_Macro_Recall", 
  "Silhouette_Sim", "PVE_Group_Sim", "pseudotime_correlation", 
  "cross_modal_accuracy", "cross_modal_F1", "match_at_1", "rv_coefficient", "module_correlation_r"
)

# ==============================================================================
# Helper Functions: Robust Matrix and Label Reading
# ==============================================================================
as_sparse_or_matrix <- function(m) {
  if (is.null(m)) return(NULL)
  if (inherits(m, "SingleCellExperiment") && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    m <- SummarizedExperiment::assay(m, "counts")
  } else if (inherits(m, "Seurat")) {
    m <- tryCatch(m[["RNA"]]$counts, error = function(e) tryCatch(Seurat::GetAssayData(m, slot = "counts"), error = function(e2) m))
  }
  if (inherits(m, "sparseMatrix")) return(m)
  if (is.data.frame(m)) {
    if (ncol(m) > 1 && (is.character(m[[1]]) || is.factor(m[[1]]))) {
      rownames(m) <- make.unique(as.character(m[[1]]))
      m <- m[, -1, drop = FALSE]
    }
    m <- data.matrix(m)
  }
  if (is.matrix(m)) {
    if (!is.numeric(m)) {
      suppressWarnings(mode(m) <- "numeric")
      m[is.na(m)] <- 0
    }
    n_tot <- as.numeric(nrow(m)) * as.numeric(ncol(m))
    if (!is.na(n_tot) && n_tot > 50000) {
      s_idx <- sample.int(n_tot, min(10000, n_tot))
      if (mean(m[s_idx] == 0, na.rm = TRUE) > 0.3) {
        return(Matrix::Matrix(m, sparse = TRUE))
      }
    }
  }
  m
}

read_uploaded_matrix <- function(file_path, file_name) {
  if (is.null(file_path) || !file.exists(file_path)) return(NULL)
  ext <- tolower(tools::file_ext(file_name))
  
  if (ext == "rds") {
    obj <- readRDS(file_path)
    if (inherits(obj, "sparseMatrix")) {
      return(obj)
    } else if (is.matrix(obj) || is.data.frame(obj)) {
      return(as_sparse_or_matrix(obj))
    } else if (inherits(obj, "SingleCellExperiment") && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      return(as_sparse_or_matrix(SingleCellExperiment::counts(obj)))
    } else if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
      return(as_sparse_or_matrix(tryCatch(obj[["RNA"]]$counts, error = function(e) Seurat::GetAssayData(obj, slot = "counts"))))
    } else if (is.list(obj) && length(obj) > 0) {
      candidate_names <- c("new_count", "counts", "sim_counts", "simulated", "sim_rna", "sim_atac", "mat")
      found <- intersect(candidate_names, names(obj))
      if (length(found) > 0) return(as_sparse_or_matrix(obj[[found[1]]]))
      for (item in obj) {
        if (is.matrix(item) || inherits(item, "Matrix") || is.data.frame(item)) {
          return(as_sparse_or_matrix(item))
        }
      }
      return(as_sparse_or_matrix(obj[[1]]))
    } else {
      stop("Unsupported RDS format. Please provide a count matrix or data.frame.")
    }
  } else if (ext == "csv") {
    df <- tryCatch(
      utils::read.csv(file_path, row.names = 1, check.names = FALSE),
      error = function(e) utils::read.csv(file_path, check.names = FALSE)
    )
    return(as_sparse_or_matrix(df))
  } else if (ext %in% c("tsv", "txt")) {
    df <- tryCatch(
      utils::read.table(file_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE),
      error = function(e) utils::read.table(file_path, sep = "\t", header = TRUE, check.names = FALSE)
    )
    return(as_sparse_or_matrix(df))
  } else if (ext %in% c("xlsx", "xls")) {
    if (requireNamespace("readxl", quietly = TRUE)) {
      df <- as.data.frame(readxl::read_excel(file_path))
      if (!is.null(df) && nrow(df) > 0 && ncol(df) > 1) {
        if (is.character(df[[1]])) {
          rownames(df) <- make.unique(as.character(df[[1]]))
          df <- df[, -1, drop = FALSE]
        }
      }
      return(as_sparse_or_matrix(df))
    } else {
      stop("readxl package required to read Excel matrix files.")
    }
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
    if (is.data.frame(obj)) {
      if (ncol(obj) >= 2) {
        c1_name <- tolower(colnames(obj)[1])
        c2_name <- tolower(colnames(obj)[2])
        if (any(grepl("cell|barcode|id", c1_name)) || any(grepl("type|cluster|batch|group|label|state", c2_name))) {
          res <- as.factor(obj[[2]])
          names(res) <- as.character(obj[[1]])
          return(res)
        }
      }
      if (ncol(obj) >= 1) return(as.factor(obj[[1]])) else return(NULL)
    }
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- if (ext == "csv") "," else "\t"
    df <- tryCatch(utils::read.table(file_path, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (!is.null(df) && is.data.frame(df) && ncol(df) >= 1) {
      if (ncol(df) >= 2) {
        c1_name <- tolower(colnames(df)[1])
        c2_name <- tolower(colnames(df)[2])
        if (any(grepl("cell|barcode|id", c1_name)) || any(grepl("type|cluster|batch|group|label|state", c2_name))) {
          res <- as.factor(df[[2]])
          names(res) <- as.character(df[[1]])
          return(res)
        }
      }
      return(as.factor(df[[1]]))
    }
  } else if (ext %in% c("xlsx", "xls")) {
    if (requireNamespace("readxl", quietly = TRUE)) {
      df <- tryCatch(as.data.frame(readxl::read_excel(file_path)), error = function(e) NULL)
      if (!is.null(df) && is.data.frame(df) && ncol(df) >= 1) {
        if (ncol(df) >= 2) {
          c1_name <- tolower(colnames(df)[1])
          c2_name <- tolower(colnames(df)[2])
          if (any(grepl("cell|barcode|id", c1_name)) || any(grepl("type|cluster|batch|group|label|state", c2_name))) {
            res <- as.factor(df[[2]])
            names(res) <- as.character(df[[1]])
            return(res)
          }
        }
        return(as.factor(df[[1]]))
      }
    }
  }
  return(NULL)
}

# ==============================================================================
# Robust Pre-Flight Input Validation & Sanitization Engine
# ==============================================================================
validate_and_sanitize_matrix <- function(mat, matrix_name = "Count Matrix", min_cells = 2, min_features = 2) {
  if (is.null(mat)) {
    return(list(valid = FALSE, error = paste0(matrix_name, ": Input is NULL or missing.")))
  }
  if (is.data.frame(mat)) {
    mat <- as_sparse_or_matrix(mat)
  }
  if (!is.matrix(mat) && !inherits(mat, "Matrix")) {
    return(list(valid = FALSE, error = paste0(matrix_name, ": Input must be a 2D matrix or data frame.")))
  }
  
  nr <- nrow(mat)
  nc <- ncol(mat)
  if (is.null(nr) || is.null(nc) || nr < min_features || nc < min_cells) {
    return(list(valid = FALSE, error = sprintf(
      "%s has insufficient dimensions (%d features x %d cells). Minimum required is %d features and %d cells.",
      matrix_name, if(is.null(nr)) 0 else nr, if(is.null(nc)) 0 else nc, min_features, min_cells
    )))
  }
  
  # Ensure rownames and colnames exist to prevent out-of-bounds on dimension lookup
  if (is.null(rownames(mat)) || any(is.na(rownames(mat))) || any(rownames(mat) == "")) {
    def_r <- paste0("Feature_", seq_len(nr))
    rownames(mat) <- make.unique(if (!is.null(rownames(mat))) ifelse(is.na(rownames(mat)) | rownames(mat) == "", def_r, rownames(mat)) else def_r)
  } else if (any(duplicated(rownames(mat)))) {
    rownames(mat) <- make.unique(rownames(mat))
  }
  
  if (is.null(colnames(mat)) || any(is.na(colnames(mat))) || any(colnames(mat) == "")) {
    def_c <- paste0("Cell_", seq_len(nc))
    colnames(mat) <- make.unique(if (!is.null(colnames(mat))) ifelse(is.na(colnames(mat)) | colnames(mat) == "", def_c, colnames(mat)) else def_c)
  } else if (any(duplicated(colnames(mat)))) {
    colnames(mat) <- make.unique(colnames(mat))
  }
  
  # Sanitize non-finite values (NA, NaN, Inf, -Inf)
  if (inherits(mat, "dgCMatrix")) {
    if (any(is.na(mat@x)) || any(is.nan(mat@x)) || any(is.infinite(mat@x))) {
      mat@x[is.na(mat@x) | is.nan(mat@x) | is.infinite(mat@x)] <- 0
    }
  } else {
    if (any(is.na(mat)) || any(is.nan(mat)) || any(is.infinite(mat))) {
      mat[is.na(mat) | is.nan(mat) | is.infinite(mat)] <- 0
    }
  }
  
  list(valid = TRUE, mat = mat)
}

safe_read_uploaded_matrix <- function(file_path, file_name, label = "Matrix") {
  if (is.null(file_path) || !file.exists(file_path)) {
    return(list(success = FALSE, error = sprintf("File for %s does not exist on the filesystem.", label)))
  }
  res <- tryCatch({
    mat <- read_uploaded_matrix(file_path, file_name)
    val <- validate_and_sanitize_matrix(mat, matrix_name = label)
    if (!val$valid) return(list(success = FALSE, error = val$error))
    list(success = TRUE, mat = val$mat)
  }, error = function(e) {
    list(success = FALSE, error = sprintf("Error reading %s ('%s'): %s", label, file_name, e$message))
  })
  res
}

safe_align_annotation <- function(ann_vec, mat, label = "Cell Types") {
  if (is.null(ann_vec) || length(ann_vec) == 0 || is.null(mat)) return(NULL)
  ann_vec <- as.character(ann_vec)
  nc <- ncol(mat)
  
  # If vector has names matching matrix colnames
  if (!is.null(names(ann_vec)) && !is.null(colnames(mat))) {
    common <- intersect(colnames(mat), names(ann_vec))
    if (length(common) >= 2) {
      aligned <- ann_vec[colnames(mat)]
      aligned[is.na(aligned)] <- "Unknown"
      names(aligned) <- colnames(mat)
      return(as.factor(aligned))
    }
  }
  
  # If length matches exactly
  if (length(ann_vec) == nc) {
    names(ann_vec) <- colnames(mat)
    return(as.factor(ann_vec))
  }
  
  warning(sprintf("[scSimEval Guard] %s length (%d) does not match cell count (%d); setting to NULL to prevent out-of-bounds error.", label, length(ann_vec), nc), call. = FALSE)
  return(NULL)
}

# Extract comprehensive dataset structural & biological properties
extract_dataset_summary <- function(mat, role = "Biological Reference", method_name = "Empirical Reference", modality = "scRNA-seq", cell_types = NULL, batch_info = NULL) {
  if (is.null(mat)) return(NULL)
  if (inherits(mat, "SingleCellExperiment") && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    mat <- SingleCellExperiment::counts(mat)
  } else if (inherits(mat, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
    mat <- Seurat::GetAssayData(mat, slot = "counts")
  }
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }
  
  n_cells <- ncol(mat)
  n_feats <- nrow(mat)
  if (is.null(n_cells) || is.null(n_feats) || n_cells == 0 || n_feats == 0) return(NULL)
  
  # Sparsity
  sparsity_pct <- if (inherits(mat, "sparseMatrix")) {
    (1 - (Matrix::nnzero(mat) / (as.numeric(n_cells) * as.numeric(n_feats)))) * 100
  } else {
    (sum(mat == 0) / (as.numeric(n_cells) * as.numeric(n_feats))) * 100
  }
  
  # Cell Types / Biological Groups
  if (!is.null(cell_types) && length(cell_types) > 0) {
    clean_ct <- stats::na.omit(as.character(cell_types))
    u_ct <- unique(clean_ct)
    n_ct <- length(u_ct)
    ct_str <- if (n_ct > 0) formatC(n_ct, format = "d", big.mark = ",") else "Not Provided"
  } else {
    ct_str <- "Not Provided"
  }
  
  # Batches / Technical Confounders
  if (!is.null(batch_info) && length(batch_info) > 0) {
    clean_b <- stats::na.omit(as.character(batch_info))
    u_b <- unique(clean_b)
    n_b <- length(u_b)
    b_str <- if (n_b > 0) formatC(n_b, format = "d", big.mark = ",") else "Not Provided"
  } else {
    b_str <- "Not Provided"
  }
  
  # Library size & detected features
  col_s <- if (inherits(mat, "Matrix")) Matrix::colSums(mat) else colSums(mat, na.rm = TRUE)
  med_lib <- stats::median(col_s, na.rm = TRUE)
  col_det <- if (inherits(mat, "dgCMatrix")) diff(mat@p) else if (inherits(mat, "sparseMatrix")) Matrix::colSums(mat > 0) else colSums(mat > 0, na.rm = TRUE)
  med_det <- stats::median(col_det, na.rm = TRUE)
  mean_expr <- if (inherits(mat, "sparseMatrix")) {
    sum(mat) / (as.numeric(n_cells) * as.numeric(n_feats))
  } else {
    mean(mat, na.rm = TRUE)
  }
  
  data.frame(
    "Dataset / Simulator" = method_name,
    "Role" = role,
    "Modality" = modality,
    "Cells (N)" = formatC(n_cells, format = "d", big.mark = ","),
    "Features (P)" = formatC(n_feats, format = "d", big.mark = ","),
    "Sparsity" = sprintf("%.2f%%", sparsity_pct),
    "Cell Types (Groups)" = ct_str,
    "Batches" = b_str,
    "Median Lib Size" = formatC(round(med_lib, 1), format = "f", digits = 1, big.mark = ","),
    "Median Detected Features" = formatC(round(med_det, 0), format = "d", big.mark = ","),
    "Mean Expression" = formatC(round(mean_expr, 3), format = "f", digits = 3),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Helper to summarize demo dataset at launch
init_demo_summary <- function(demo_obj) {
  if (is.null(demo_obj)) return(NULL)
  if (!is.null(demo_obj$dataset_summary_df)) return(demo_obj$dataset_summary_df)
  NULL
}

# Helper to initialize demo simulated matrices list
init_demo_sim_matrices <- function(demo_obj) {
  if (is.null(demo_obj)) return(list())
  res <- list()
  if (!is.null(demo_obj$example_rice_bud) && !is.null(demo_obj$example_rice_bud$simulated)) {
    for (m in names(demo_obj$example_rice_bud$simulated)) {
      res[[paste0("Rice Bud - ", m, " (RNA)")]] <- demo_obj$example_rice_bud$simulated[[m]]$rna
      res[[paste0("Rice Bud - ", m, " (ATAC)")]] <- demo_obj$example_rice_bud$simulated[[m]]$atac
    }
  }
  if (!is.null(demo_obj$example_soybean_root) && !is.null(demo_obj$example_soybean_root$simulated)) {
    for (m in names(demo_obj$example_soybean_root$simulated)) {
      res[[paste0("Soybean Root - ", m, " (RNA)")]] <- demo_obj$example_soybean_root$simulated[[m]]$rna
      res[[paste0("Soybean Root - ", m, " (ATAC)")]] <- demo_obj$example_soybean_root$simulated[[m]]$atac
    }
  }
  res
}

# ------------------------------------------------------------------------------
# Cell Embeddings (UMAP, t-SNE, PCA) & Quality Metrics Helpers
# ------------------------------------------------------------------------------
compute_dataset_embeddings <- function(reference,
                                       simulated,
                                       reduction = c("umap", "tsne", "pca"),
                                       n_pcs = 30,
                                       perplexity = 30,
                                       n_neighbors = 15,
                                       min_dist = 0.3,
                                       seed = 42,
                                       cell_types = NULL,
                                       batch = NULL) {
  reduction <- match.arg(reduction)
  set.seed(seed)
  
  extract_counts <- function(obj) {
    if (is.null(obj)) return(NULL)
    if (inherits(obj, "SingleCellExperiment") && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      SingleCellExperiment::counts(obj)
    } else if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
      Seurat::GetAssayData(obj, slot = "counts")
    } else {
      obj
    }
  }
  
  ref_mat <- extract_counts(reference)
  if (is.null(ref_mat) || length(dim(ref_mat)) < 2) {
    stop("A valid 2D reference count matrix must be provided.", call. = FALSE)
  }
  
  sim_list <- if (is.list(simulated) && !is.data.frame(simulated) && !inherits(simulated, "dgCMatrix")) {
    simulated
  } else {
    list("Simulated" = simulated)
  }
  sim_list <- lapply(sim_list, extract_counts)
  
  if (is.null(names(sim_list)) || any(names(sim_list) == "")) {
    names(sim_list) <- paste0("Simulator_", seq_along(sim_list))
  }
  
  all_datasets <- c(list("Reference" = ref_mat), sim_list)
  dataset_names <- names(all_datasets)
  results <- list()
  
  # Extract biological cell types for the empirical reference
  ref_ct_actual <- NULL
  if (!is.null(cell_types)) {
    if (is.list(cell_types) && !is.data.frame(cell_types) && !is.factor(cell_types)) {
      if (!is.null(cell_types[["Reference"]]) && length(cell_types[["Reference"]]) == ncol(ref_mat)) {
        ref_ct_actual <- as.character(cell_types[["Reference"]])
      } else {
        for (item in cell_types) {
          if (length(item) == ncol(ref_mat)) {
            ref_ct_actual <- as.character(item)
            break
          }
        }
      }
    } else if (is.data.frame(cell_types)) {
      v <- if (!is.null(cell_types$CellType)) cell_types$CellType else if (!is.null(cell_types$Cell_Type)) cell_types$Cell_Type else cell_types[[1]]
      if (length(v) == ncol(ref_mat)) ref_ct_actual <- as.character(v)
    } else if (length(cell_types) == ncol(ref_mat)) {
      ref_ct_actual <- as.character(cell_types)
    }
  }
  ref_pca_centroids <- list()
  
  for (dname in dataset_names) {
    mat <- all_datasets[[dname]]
    if (is.null(mat) || length(dim(mat)) < 2) next
    
    n_cells <- ncol(mat)
    n_feats <- nrow(mat)
    
    if (is.null(n_cells) || is.null(n_feats) || is.na(n_cells) || is.na(n_feats) || n_cells < 3 || n_feats < 3) next
    
    # 1. Total library size and detected features
    libs <- if (inherits(mat, "sparseMatrix")) Matrix::colSums(mat) else colSums(mat)
    det_feats <- if (inherits(mat, "dgCMatrix")) diff(mat@p) else if (inherits(mat, "sparseMatrix")) Matrix::colSums(mat > 0) else colSums(mat > 0)
    
    # 2. Fast Variance-based feature selection (prioritize top 2,000 HVGs BEFORE full densification)
    top_n <- min(2000, n_feats)
    if (n_feats > top_n) {
      feat_vars <- if (exists("fast_row_vars", mode = "function")) {
        fast_row_vars(mat)
      } else if (inherits(mat, "dgCMatrix")) {
        n_c <- ncol(mat)
        rm <- Matrix::rowMeans(mat)
        m2 <- mat; m2@x <- m2@x^2
        v <- (Matrix::rowMeans(m2) - rm^2) * (n_c / (n_c - 1))
        pmax(0, as.numeric(v))
      } else {
        apply(mat, 1, stats::var)
      }
      feat_vars[is.na(feat_vars)] <- 0
      top_idx <- order(feat_vars, decreasing = TRUE)[seq_len(top_n)]
      mat_sub <- mat[top_idx, , drop = FALSE]
    } else {
      mat_sub <- mat
    }
    
    # 3. Library size scaling and log-transformation on top features only (drastically reduces RAM and runtime)
    scale_factor <- stats::median(libs[libs > 0])
    if (is.na(scale_factor) || scale_factor <= 0) scale_factor <- 10000
    
    sub_dense <- as.matrix(mat_sub)
    sub_mat <- log1p(sweep(sub_dense, 2, libs / scale_factor, "/"))
    sub_mat[is.na(sub_mat) | is.infinite(sub_mat)] <- 0
    
    # 4. Principal Component Analysis (PCA)
    k_pc <- min(n_pcs, n_cells - 1, nrow(sub_mat) - 1)
    if (k_pc < 2) k_pc <- 2
    
    pca_res <- if (exists("fast_pca", mode = "function")) {
      list(x = fast_pca(t(sub_mat), n_pcs = k_pc, scale = FALSE))
    } else if (requireNamespace("irlba", quietly = TRUE) && k_pc < (n_cells - 2) && k_pc < 0.5 * min(n_cells, nrow(sub_mat))) {
      tryCatch(
        suppressWarnings(irlba::prcomp_irlba(t(sub_mat), n = k_pc, center = TRUE, scale. = FALSE)),
        error = function(e) stats::prcomp(t(sub_mat), center = TRUE, scale. = FALSE)
      )
    } else {
      stats::prcomp(t(sub_mat), center = TRUE, scale. = FALSE)
    }
    pca_coords <- pca_res$x[, seq_len(min(k_pc, ncol(pca_res$x))), drop = FALSE]
    
    # 5. Non-linear Dimension Reduction (UMAP / t-SNE / PCA)
    dim1 <- if (ncol(pca_coords) >= 1) pca_coords[, 1] else rep(0, n_cells)
    dim2 <- if (ncol(pca_coords) >= 2) pca_coords[, 2] else dim1
    
    if (reduction == "tsne") {
      perp <- min(perplexity, floor((n_cells - 1) / 3))
      if (perp < 2) perp <- 2
      if (requireNamespace("Rtsne", quietly = TRUE)) {
        tryCatch({
          tsne_out <- Rtsne::Rtsne(pca_coords, perplexity = perp, check_duplicates = FALSE, pca = FALSE)
          if (!is.null(tsne_out$Y) && ncol(tsne_out$Y) >= 2) {
            dim1 <- tsne_out$Y[, 1]
            dim2 <- tsne_out$Y[, 2]
          }
        }, error = function(e) {
          warning("t-SNE computation failed: ", e$message, "; falling back to PCA coordinates.", call. = FALSE)
        })
      }
    } else if (reduction == "umap") {
      n_neigh <- min(n_neighbors, n_cells - 1)
      if (n_neigh < 2) n_neigh <- 2
      if (requireNamespace("uwot", quietly = TRUE)) {
        tryCatch({
          umap_out <- uwot::umap(pca_coords, n_neighbors = n_neigh, min_dist = min_dist, seed = seed)
          if (!is.null(umap_out) && ncol(umap_out) >= 2) {
            dim1 <- umap_out[, 1]
            dim2 <- umap_out[, 2]
          }
        }, error = function(e) {
          warning("UMAP computation failed: ", e$message, "; falling back to PCA coordinates.", call. = FALSE)
        })
      }
    }
    
    # 6. Unsupervised Cluster Recovery (k-means on PCA)
    k_clust <- if (!is.null(ref_ct_actual)) {
      length(unique(stats::na.omit(ref_ct_actual)))
    } else if (!is.null(cell_types) && !is.list(cell_types)) {
      length(unique(stats::na.omit(as.character(cell_types))))
    } else {
      3
    }
    k_clust <- max(2, min(k_clust, n_cells - 1))
    km_fit <- tryCatch(
      stats::kmeans(pca_coords, centers = k_clust, nstart = 5),
      error = function(e) list(cluster = rep(1, n_cells), centers = matrix(0, 1, ncol(pca_coords)))
    )
    
    # 7. Exact Biological Cell Type Alignment & Nearest-Neighbor Transfer
    ct_vec <- NULL
    
    # 1. Direct list lookup
    if (is.list(cell_types) && !is.data.frame(cell_types) && !is.factor(cell_types)) {
      if (!is.null(cell_types[[dname]]) && length(cell_types[[dname]]) == n_cells) {
        ct_vec <- as.character(cell_types[[dname]])
      } else {
        m_idx <- match(tolower(dname), tolower(names(cell_types)))
        if (!is.na(m_idx) && length(cell_types[[m_idx]]) == n_cells) {
          ct_vec <- as.character(cell_types[[m_idx]])
        }
      }
    }
    
    # 2. Match by cell barcode/name
    if (is.null(ct_vec) && !is.null(cell_types) && !is.null(names(cell_types)) && !is.null(colnames(mat))) {
      m_cells <- match(colnames(mat), names(cell_types))
      if (sum(!is.na(m_cells)) >= 0.7 * n_cells) {
        ct_vec <- as.character(cell_types[m_cells])
        ct_vec[is.na(ct_vec)] <- "Unknown"
      }
    }
    
    # 3. Exact cell count match
    if (is.null(ct_vec) && !is.null(cell_types) && !is.list(cell_types) && length(cell_types) == n_cells) {
      ct_vec <- as.character(cell_types)
    }
    
    # 4. Reference dataset fallback to ref_ct_actual
    if (is.null(ct_vec) && dname == "Reference" && !is.null(ref_ct_actual) && length(ref_ct_actual) == n_cells) {
      ct_vec <- ref_ct_actual
    }
    
    # 5. If this is Reference, compute reference PCA centroids for biological cell types
    if (dname == "Reference" && !is.null(ct_vec) && !all(grepl("^Cluster ", ct_vec))) {
      ref_pca_centroids <- list()
      for (u_ct in unique(ct_vec)) {
        idx_u <- which(ct_vec == u_ct)
        if (length(idx_u) > 0) {
          ref_pca_centroids[[u_ct]] <- colMeans(pca_coords[idx_u, , drop = FALSE])
        }
      }
    }
    
    # 6. Simulator dataset: Transfer biological cell type labels from Reference via PCA centroid proximity
    if (is.null(ct_vec) && dname != "Reference" && length(ref_pca_centroids) >= 2) {
      ct_vec <- tryCatch({
        km_centers <- km_fit$centers
        u_ref_types <- names(ref_pca_centroids)
        clust_to_type <- character(nrow(km_centers))
        for (ci in seq_len(nrow(km_centers))) {
          c_vec <- km_centers[ci, , drop = TRUE]
          dists <- vapply(u_ref_types, function(rt) {
            cent <- ref_pca_centroids[[rt]]
            k_len <- min(length(c_vec), length(cent))
            sqrt(sum((c_vec[seq_len(k_len)] - cent[seq_len(k_len)])^2))
          }, numeric(1))
          clust_to_type[ci] <- u_ref_types[which.min(dists)]
        }
        clust_to_type[km_fit$cluster]
      }, error = function(e) NULL)
    }
    
    # 7. Final fallback: Use ref_ct_actual if length matches, else "Cluster X"
    if (is.null(ct_vec) || any(is.na(ct_vec)) || length(ct_vec) != n_cells) {
      if (!is.null(ref_ct_actual) && length(ref_ct_actual) == n_cells) {
        ct_vec <- ref_ct_actual
      } else {
        ct_vec <- paste0("Cluster ", km_fit$cluster)
      }
    }
    
    b_vec <- if (!is.null(batch)) {
      if (is.list(batch) && !is.null(batch[[dname]]) && length(batch[[dname]]) == n_cells) {
        as.character(batch[[dname]])
      } else if (length(batch) == n_cells) {
        as.character(batch)
      } else {
        "Batch 1"
      }
    } else {
      "Batch 1"
    }
    
    c_names <- if (!is.null(colnames(mat))) colnames(mat) else paste0("Cell_", seq_len(n_cells))
    
    df <- data.frame(
      Cell_ID = c_names,
      Dim1 = as.numeric(dim1),
      Dim2 = as.numeric(dim2),
      Dataset = dname,
      Role = ifelse(dname == "Reference", "Reference", "Simulated"),
      Dataset_Type = ifelse(dname == "Reference", "Reference", "Simulated"),
      Cell_Type = as.character(ct_vec),
      Cluster = factor(paste0("Cluster_", km_fit$cluster)),
      Library_Size = as.numeric(libs),
      Detected_Features = as.numeric(det_feats),
      Batch = as.character(b_vec),
      stringsAsFactors = FALSE
    )
    results[[dname]] <- df
  }
  
  if (length(results) == 0) return(data.frame())
  
  combined_df <- do.call(rbind, results)
  combined_df$Dataset <- factor(combined_df$Dataset, levels = dataset_names)
  rownames(combined_df) <- NULL
  combined_df
}

plot_dataset_embeddings <- function(embedding_data,
                                    reduction = c("umap", "tsne", "pca"),
                                    layout = c("facet", "side_by_side", "overlay"),
                                    color_by = c("cell_type", "cluster", "library_size", "dataset", "batch"),
                                    selected_methods = NULL,
                                    pt_size = 0.8,
                                    alpha = 0.75,
                                    palette = "Set1",
                                    title = NULL,
                                    base_size = 12) {
  reduction <- match.arg(reduction)
  layout <- match.arg(layout)
  color_by <- match.arg(color_by)
  
  if (is.null(embedding_data) || nrow(embedding_data) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No embedding coordinates available to plot.", size = 5, color = "#64748B") +
        ggplot2::theme_void()
    )
  }
  
  df <- embedding_data
  red_label <- switch(reduction, "umap" = "UMAP", "tsne" = "t-SNE", "pca" = "PCA")
  
  # Filter selected methods if requested
  if (!is.null(selected_methods) && length(selected_methods) > 0) {
    keep_datasets <- c("Reference", selected_methods)
    df_datasets <- as.character(df$Dataset)
    base_datasets <- sub("(_Dataset.*| \\(Dataset.*\\))", "", df_datasets)
    keep_idx <- df_datasets %in% keep_datasets | base_datasets %in% selected_methods
    df <- df[keep_idx, , drop = FALSE]
  }
  
  if (layout == "side_by_side") {
    sim_levels <- setdiff(unique(as.character(df$Dataset)), "Reference")
    chosen_sim <- if (!is.null(selected_methods) && length(selected_methods) == 1 &&
                      any(selected_methods == sim_levels | selected_methods == sub("(_Dataset.*| \\(Dataset.*\\))", "", sim_levels))) {
      sim_levels[which(selected_methods == sim_levels | selected_methods == sub("(_Dataset.*| \\(Dataset.*\\))", "", sim_levels))[1]]
    } else if (length(sim_levels) > 0) {
      sim_levels[1]
    } else {
      NULL
    }
    if (!is.null(chosen_sim)) {
      df <- df[df$Dataset %in% c("Reference", chosen_sim), , drop = FALSE]
    }
  }
  
  df$Facet_Label <- ifelse(
    df$Dataset == "Reference" | df$Role == "Reference",
    "Reference",
    paste0(df$Dataset, " (", df$Role, ")")
  )
  df$Facet_Label <- factor(df$Facet_Label, levels = unique(df$Facet_Label))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2))
  
  if (color_by == "cell_type") {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Cell_Type), size = pt_size, alpha = alpha) +
      ggplot2::labs(color = "Cell Type")
    n_colors <- length(unique(stats::na.omit(df$Cell_Type)))
    if (n_colors <= 9 && requireNamespace("RColorBrewer", quietly = TRUE)) {
      p <- p + ggplot2::scale_color_brewer(palette = palette)
    } else {
      cols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(8, max(3, n_colors)), palette))(n_colors)
      p <- p + ggplot2::scale_color_manual(values = cols)
    }
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3.5, alpha = 1)))
  } else if (color_by == "cluster") {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Cluster), size = pt_size, alpha = alpha) +
      ggplot2::labs(color = "Cluster")
    n_colors <- length(unique(stats::na.omit(df$Cluster)))
    if (n_colors <= 9 && requireNamespace("RColorBrewer", quietly = TRUE)) {
      p <- p + ggplot2::scale_color_brewer(palette = "Set2")
    } else {
      cols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(8, max(3, n_colors)), "Set2"))(n_colors)
      p <- p + ggplot2::scale_color_manual(values = cols)
    }
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3.5, alpha = 1)))
  } else if (color_by == "library_size") {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Library_Size), size = pt_size, alpha = alpha) +
      ggplot2::scale_color_viridis_c(option = "plasma", name = "Library Size")
  } else if (color_by == "batch") {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Batch), size = pt_size, alpha = alpha) +
      ggplot2::labs(color = "Batch")
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Dataset), size = pt_size, alpha = alpha) +
      ggplot2::labs(color = "Dataset")
  }
  
  p <- p +
    ggplot2::labs(
      title = NULL,
      subtitle = NULL,
      x = paste(red_label, "1"),
      y = paste(red_label, "2")
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = base_size, color = "#1E293B"),
      axis.text = ggplot2::element_text(size = base_size * 0.85, color = "#475569"),
      strip.text = ggplot2::element_text(face = "bold", size = base_size * 0.95, color = "#0F172A"),
      strip.background = ggplot2::element_rect(fill = "#F8FAFC", color = "#CBD5E1", linewidth = 0.5),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold", size = base_size * 0.90, color = "#1E293B"),
      legend.text = ggplot2::element_text(size = base_size * 0.85, color = "#334155"),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "#CBD5E1", fill = NA, linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#F1F5F9", linewidth = 0.5)
    )
  
  if (layout == "facet") {
    p <- p + ggplot2::facet_wrap(~ Facet_Label, scales = "fixed")
  } else if (layout == "side_by_side") {
    p <- p + ggplot2::facet_wrap(~ Facet_Label, nrow = 1, scales = "fixed")
  }
  
  p
}

compute_embedding_quality_metrics <- function(embedding_data, metric_coords = c("Dim1", "Dim2")) {
  if (is.null(embedding_data) || nrow(embedding_data) == 0) return(data.frame())
  
  datasets <- unique(as.character(embedding_data$Dataset))
  rows <- list()
  
  for (d in datasets) {
    sub <- embedding_data[embedding_data$Dataset == d, , drop = FALSE]
    n_c <- nrow(sub)
    
    # 1. Silhouette score on cell types
    sil_val <- NA_real_
    if (length(unique(stats::na.omit(sub$Cell_Type))) > 1 && n_c >= 4) {
      coords <- as.matrix(sub[, metric_coords, drop = FALSE])
      labels <- as.integer(factor(sub$Cell_Type))
      idx_s <- sample(seq_len(n_c), size = min(2000, n_c))
      dmat <- stats::dist(coords[idx_s, , drop = FALSE])
      if (requireNamespace("cluster", quietly = TRUE)) {
        sil_obj <- cluster::silhouette(labels[idx_s], dmat)
        sil_val <- mean(sil_obj[, 3], na.rm = TRUE)
      }
    }
    
    # 2. Adjusted Rand Index (ARI) of cluster recovery vs ground truth
    ari_val <- NA_real_
    if (!is.null(sub$Cluster) && length(unique(stats::na.omit(sub$Cluster))) > 1 &&
        length(unique(stats::na.omit(sub$Cell_Type))) > 1) {
      if (requireNamespace("mclust", quietly = TRUE)) {
        ari_val <- mclust::adjustedRandIndex(sub$Cluster, sub$Cell_Type)
      }
    }
    
    # 3. Distribution metrics
    m_umi <- mean(sub$Library_Size, na.rm = TRUE)
    m_det <- mean(sub$Detected_Features, na.rm = TRUE)
    role_val <- if (!is.null(sub$Role)) unique(sub$Role)[1] else unique(sub$Dataset_Type)[1]
    
    rows[[d]] <- data.frame(
      Dataset = d,
      Role = role_val,
      "Cells (N)" = n_c,
      "Mean Silhouette" = round(sil_val, 4),
      "ARI (Cluster Fidelity)" = round(ari_val, 4),
      "Mean Library Size" = round(m_umi, 1),
      "Mean Detected Features" = round(m_det, 1),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  res_df <- do.call(rbind, rows)
  rownames(res_df) <- NULL
  
  # 4. Difference % relative to Reference
  ref_row <- res_df[res_df$Dataset == "Reference", , drop = FALSE]
  if (nrow(ref_row) > 0) {
    ref_umi <- ref_row[["Mean Library Size"]][1]
    res_df[["Library Size Diff (%)"]] <- round(100 * abs(res_df[["Mean Library Size"]] - ref_umi) / max(ref_umi, 1), 2)
  }
  
  res_df
}

# Export multi-sheet Excel workbook
export_excel_workbook <- function(file, benchmark_df, leaderboard_df = NULL, dataset_summary_df = NULL) {
  sheets <- list(All_Benchmark_Metrics = benchmark_df)
  if (!is.null(dataset_summary_df) && nrow(dataset_summary_df) > 0) {
    sheets$Dataset_Properties <- dataset_summary_df
  }
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

# Export High-Res Publication JPEG (600 DPI)
export_single_jpeg <- function(file, plot_obj, width = 14, height = 9, dpi = 600) {
  ggplot2::ggsave(file, plot = plot_obj, device = "jpeg", width = width, height = height, dpi = dpi)
}

# Export Vectorized Publication PDF
export_single_pdf <- function(file, plot_obj, width = 14, height = 9) {
  grDevices::pdf(file, width = width, height = height)
  try(print(plot_obj), silent = TRUE)
  grDevices::dev.off()
}

# Generate Multi-Page PDF Report
generate_all_plots_pdf <- function(file, benchmark_data, toy_ref = NULL, toy_sim = NULL, sim_matrices = NULL, cell_types = NULL, batch = NULL) {
  grDevices::pdf(file, width = 14, height = 9, onefile = TRUE)
  try(print(plot_benchmark_bubble_matrix(benchmark_data, base_size = 9.5, show_missing_dots = FALSE)), silent = TRUE)
  try(print(plot_evaluation_summary(benchmark_data, base_size = 12)), silent = TRUE)
  try(print(plot_scalability_benchmark(benchmark_data, base_size = 12)), silent = TRUE)
  try(print(plot_metric_boxplots(benchmark_data, base_size = 11)), silent = TRUE)
  try(print(plot_metric_heatmap(benchmark_data, base_size = 11)), silent = TRUE)
  try(print(plot_metric_pca(benchmark_data, base_size = 12)), silent = TRUE)
  try(print(plot_metric_mds(benchmark_data, base_size = 12)), silent = TRUE)
  if (!is.null(toy_ref) && !is.null(toy_sim)) {
    try(print(plot_distribution_qc(toy_ref, toy_sim, base_size = 11)), silent = TRUE)
    s_list <- if (!is.null(sim_matrices) && length(sim_matrices) > 0) sim_matrices else list("Simulated" = toy_sim)
    for (red in c("umap", "tsne", "pca")) {
      try({
        emb <- compute_dataset_embeddings(toy_ref, s_list, reduction = red, cell_types = cell_types, batch = batch)
        print(plot_dataset_embeddings(emb, reduction = red, layout = "facet", base_size = 11))
      }, silent = TRUE)
    }
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
    tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"),
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
      
      /* Horizontal & Vertical scroll container for big bubble plot */
      .bubble-scroll-container {
        overflow-x: auto;
        overflow-y: auto;
        max-height: 720px;
        border: 1px solid #E2E8F0;
        border-radius: 8px;
        background: #FFFFFF;
        padding: 12px;
        box-shadow: inset 0 1px 3px rgba(0,0,0,0.03);
      }
      
      /* Clean scientific inputs */
      .sim-input-card {
        background: #F8FAFC;
        border-left: 4px solid #1B4F72;
        border-radius: 6px;
        padding: 12px;
        margin-bottom: 12px;
      }

      /* Team & Contact Info Cards */
      .info-card {
        background: #FFFFFF;
        border: 1px solid #CBD5E1;
        border-radius: 12px;
        padding: 18px 14px;
        margin: 12px auto;
        max-width: 290px;
        transition: all 0.25s ease;
        box-shadow: 0 2px 6px rgba(0,0,0,0.05);
        text-align: center;
      }
      .info-card:hover {
        transform: translateY(-4px);
        box-shadow: 0 8px 20px rgba(0,0,0,0.12);
        border-color: #3B82F6;
        background: #FFFFFF;
      }

      /* ===== GLOBAL APP-WIDE LOADING INDICATOR ===== */
      .app-global-loader {
        position: fixed;
        top: 14px;
        right: 22px;
        background: linear-gradient(135deg, #1E3A8A 0%, #0F172A 100%);
        color: #FFFFFF;
        padding: 8px 18px;
        border-radius: 24px;
        font-size: 0.82rem;
        font-weight: 700;
        letter-spacing: 0.03em;
        box-shadow: 0 4px 16px rgba(15, 23, 42, 0.35);
        z-index: 99999;
        display: none;
        align-items: center;
        gap: 10px;
        pointer-events: none;
        border: 1px solid rgba(255, 255, 255, 0.2);
        animation: pulseLoader 1.5s infinite alternate;
      }
      @keyframes pulseLoader {
        from { transform: scale(1); }
        to { transform: scale(1.03); }
      }
      html.shiny-busy .app-global-loader {
        display: inline-flex !important;
      }
      .shiny-plot-output.recalculating,
      .shiny-html-output.recalculating {
        opacity: 0.35 !important;
        filter: blur(1px);
        transition: all 0.2s ease;
      }
      /* ===== GENERAL DAYLIGHT TYPOGRAPHY & CONTRAST (DAYLIGHT MODE ONLY) ===== */
      html:not([data-bs-theme='dark']) body {
        color: #0F172A;
      }
      html:not([data-bs-theme='dark']) .text-muted,
      html:not([data-bs-theme='dark']) small.text-muted {
        color: #475569 !important;
      }
      
      /* ===== HIGH-CONTRAST ALERT BANNERS IN DAYLIGHT MODE ONLY ===== */
      /* Override Bootswatch Flatly default white text on light alert cards */
      html:not([data-bs-theme='dark']) .alert {
        color: #0F172A !important;
      }
      html:not([data-bs-theme='dark']) .alert-info {
        background-color: #F0F9FF !important;
        border-color: #BAE6FD !important;
        color: #0369A1 !important;
      }
      html:not([data-bs-theme='dark']) .alert-info,
      html:not([data-bs-theme='dark']) .alert-info p,
      html:not([data-bs-theme='dark']) .alert-info div,
      html:not([data-bs-theme='dark']) .alert-info span:not(.badge) {
        color: #0369A1 !important;
      }
      html:not([data-bs-theme='dark']) .alert-info b,
      html:not([data-bs-theme='dark']) .alert-info strong,
      html:not([data-bs-theme='dark']) .alert-info .alert-heading {
        color: #0C4A6E !important;
      }
      html:not([data-bs-theme='dark']) .alert-info i,
      html:not([data-bs-theme='dark']) .alert-info .fa,
      html:not([data-bs-theme='dark']) .alert-info .fas {
        color: #0284C7 !important;
      }
      html:not([data-bs-theme='dark']) .alert-info a,
      html:not([data-bs-theme='dark']) .alert-info a:visited,
      html:not([data-bs-theme='dark']) .alert-info .action-button,
      html:not([data-bs-theme='dark']) .alert-info actionLink {
        color: #0284C7 !important;
        font-weight: 700 !important;
      }

      html:not([data-bs-theme='dark']) .alert-primary {
        background-color: #EEF2FF !important;
        border-color: #C7D2FE !important;
        color: #3730A3 !important;
      }
      html:not([data-bs-theme='dark']) .alert-primary,
      html:not([data-bs-theme='dark']) .alert-primary p,
      html:not([data-bs-theme='dark']) .alert-primary div,
      html:not([data-bs-theme='dark']) .alert-primary span:not(.badge) {
        color: #3730A3 !important;
      }
      html:not([data-bs-theme='dark']) .alert-primary b,
      html:not([data-bs-theme='dark']) .alert-primary strong,
      html:not([data-bs-theme='dark']) .alert-primary .alert-heading {
        color: #1E1B4B !important;
      }
      html:not([data-bs-theme='dark']) .alert-primary i,
      html:not([data-bs-theme='dark']) .alert-primary .fa,
      html:not([data-bs-theme='dark']) .alert-primary .fas {
        color: #4F46E5 !important;
      }
      html:not([data-bs-theme='dark']) .alert-primary a,
      html:not([data-bs-theme='dark']) .alert-primary a:visited,
      html:not([data-bs-theme='dark']) .alert-primary .action-button,
      html:not([data-bs-theme='dark']) .alert-primary actionLink {
        color: #4F46E5 !important;
        font-weight: 700 !important;
      }

      html:not([data-bs-theme='dark']) .alert-indigo {
        background-color: #EEF2FF !important;
        border-color: #C7D2FE !important;
        color: #3730A3 !important;
      }
      html:not([data-bs-theme='dark']) .alert-indigo,
      html:not([data-bs-theme='dark']) .alert-indigo p,
      html:not([data-bs-theme='dark']) .alert-indigo div,
      html:not([data-bs-theme='dark']) .alert-indigo span:not(.badge) {
        color: #3730A3 !important;
      }
      html:not([data-bs-theme='dark']) .alert-indigo b,
      html:not([data-bs-theme='dark']) .alert-indigo strong {
        color: #1E1B4B !important;
      }
      html:not([data-bs-theme='dark']) .alert-indigo i,
      html:not([data-bs-theme='dark']) .alert-indigo .fa,
      html:not([data-bs-theme='dark']) .alert-indigo .fas {
        color: #4F46E5 !important;
      }

      html:not([data-bs-theme='dark']) .alert-secondary {
        background-color: #F8FAFC !important;
        border-color: #CBD5E1 !important;
        color: #334155 !important;
      }
      html:not([data-bs-theme='dark']) .alert-secondary,
      html:not([data-bs-theme='dark']) .alert-secondary p,
      html:not([data-bs-theme='dark']) .alert-secondary div,
      html:not([data-bs-theme='dark']) .alert-secondary span:not(.badge) {
        color: #334155 !important;
      }
      html:not([data-bs-theme='dark']) .alert-secondary b,
      html:not([data-bs-theme='dark']) .alert-secondary strong {
        color: #0F172A !important;
      }

      html:not([data-bs-theme='dark']) .alert-success {
        background-color: #F0FDF4 !important;
        border-color: #BBF7D0 !important;
        color: #166534 !important;
      }
      html:not([data-bs-theme='dark']) .alert-success,
      html:not([data-bs-theme='dark']) .alert-success p,
      html:not([data-bs-theme='dark']) .alert-success div,
      html:not([data-bs-theme='dark']) .alert-success span:not(.badge) {
        color: #166534 !important;
      }
      html:not([data-bs-theme='dark']) .alert-success b,
      html:not([data-bs-theme='dark']) .alert-success strong {
        color: #14532D !important;
      }

      html:not([data-bs-theme='dark']) .alert-warning {
        background-color: #FFFBEB !important;
        border-color: #FDE68A !important;
        color: #92400E !important;
      }
      html:not([data-bs-theme='dark']) .alert-warning,
      html:not([data-bs-theme='dark']) .alert-warning p,
      html:not([data-bs-theme='dark']) .alert-warning div,
      html:not([data-bs-theme='dark']) .alert-warning span:not(.badge) {
        color: #92400E !important;
      }
      html:not([data-bs-theme='dark']) .alert-warning b,
      html:not([data-bs-theme='dark']) .alert-warning strong {
        color: #78350F !important;
      }

      html:not([data-bs-theme='dark']) .alert-danger {
        background-color: #FEF2F2 !important;
        border-color: #FECACA !important;
        color: #991B1B !important;
      }
      html:not([data-bs-theme='dark']) .alert-danger,
      html:not([data-bs-theme='dark']) .alert-danger p,
      html:not([data-bs-theme='dark']) .alert-danger div,
      html:not([data-bs-theme='dark']) .alert-danger span:not(.badge) {
        color: #991B1B !important;
      }
      html:not([data-bs-theme='dark']) .alert-danger b,
      html:not([data-bs-theme='dark']) .alert-danger strong {
        color: #7F1D1D !important;
      }

      /* Alert Badges High Contrast in Daylight */
      html:not([data-bs-theme='dark']) .alert .badge.bg-primary-subtle,
      html:not([data-bs-theme='dark']) .badge.bg-primary-subtle {
        background-color: #E0E7FF !important;
        color: #3730A3 !important;
        border: 1px solid #C7D2FE !important;
        font-weight: 700 !important;
      }
      html:not([data-bs-theme='dark']) .alert .badge.bg-warning-subtle,
      html:not([data-bs-theme='dark']) .badge.bg-warning-subtle {
        background-color: #FEF3C7 !important;
        color: #92400E !important;
        border: 1px solid #FDE68A !important;
        font-weight: 700 !important;
      }
      
      /* Daylight Mode: Clean scientific confirmation & placeholder card */
      html:not([data-bs-theme='dark']) .plot-placeholder-card {
        text-align: center;
        padding: 48px 24px;
        background: #F8FAFC !important;
        border: 2px dashed #CBD5E1 !important;
        border-radius: 12px;
        margin: 24px auto;
        max-width: 680px;
        box-shadow: 0 1px 4px rgba(0,0,0,0.04);
      }
      html:not([data-bs-theme='dark']) .plot-placeholder-card h5,
      html:not([data-bs-theme='dark']) .plot-placeholder-card .placeholder-title {
        font-weight: 700;
        color: #0F172A !important;
        margin-bottom: 8px;
        font-size: 1.15rem;
      }
      html:not([data-bs-theme='dark']) .plot-placeholder-card p,
      html:not([data-bs-theme='dark']) .plot-placeholder-card .placeholder-desc {
        color: #475569 !important;
        font-size: 0.88rem !important;
        line-height: 1.55;
        margin-bottom: 20px;
      }
      
      /* ==================================================================== */
      /* COMPREHENSIVE HIGH-CONTRAST NIGHT MODE / DARK THEME STYLING          */
      /* ==================================================================== */
      
      /* 1. Global Dark Theme CSS Variables */
      [data-bs-theme='dark'] {
        --bs-body-bg: #0B1329 !important;
        --bs-body-color: #F8FAFC !important;
        --bs-emphasis-color: #FFFFFF !important;
        --bs-secondary-color: #CBD5E1 !important;
        --bs-tertiary-color: #94A3B8 !important;
        --bs-border-color: #334155 !important;
        --bs-card-bg: #0F172A !important;
        --bs-card-color: #F8FAFC !important;
        --bs-card-border-color: #334155 !important;
        --bs-card-cap-bg: #1E293B !important;
      }
      
      [data-bs-theme='dark'],
      [data-bs-theme='dark'] body {
        background-color: #0B1329 !important;
        color: #F8FAFC !important;
      }

      /* 2. Universal Light Typography in Night Mode */
      [data-bs-theme='dark'] p,
      [data-bs-theme='dark'] span:not(.badge),
      [data-bs-theme='dark'] div:not(.badge):not(.btn):not(.selectize-dropdown-content),
      [data-bs-theme='dark'] label,
      [data-bs-theme='dark'] .form-label,
      [data-bs-theme='dark'] .control-label,
      [data-bs-theme='dark'] li,
      [data-bs-theme='dark'] td,
      [data-bs-theme='dark'] th,
      [data-bs-theme='dark'] small,
      [data-bs-theme='dark'] em,
      [data-bs-theme='dark'] summary,
      [data-bs-theme='dark'] caption,
      [data-bs-theme='dark'] .help-block,
      [data-bs-theme='dark'] .form-text,
      [data-bs-theme='dark'] .shiny-input-container label {
        color: #F8FAFC !important;
      }

      /* Bold and Strong text: Pure High-Contrast White */
      [data-bs-theme='dark'] strong,
      [data-bs-theme='dark'] b,
      [data-bs-theme='dark'] .fw-bold,
      [data-bs-theme='dark'] .fw-semibold {
        color: #FFFFFF !important;
      }

      /* Headings: Crisp Pure White */
      [data-bs-theme='dark'] h1, [data-bs-theme='dark'] .h1,
      [data-bs-theme='dark'] h2, [data-bs-theme='dark'] .h2,
      [data-bs-theme='dark'] h3, [data-bs-theme='dark'] .h3,
      [data-bs-theme='dark'] h4, [data-bs-theme='dark'] .h4,
      [data-bs-theme='dark'] h5, [data-bs-theme='dark'] .h5,
      [data-bs-theme='dark'] h6, [data-bs-theme='dark'] .h6 {
        color: #FFFFFF !important;
      }

      /* Secondary / Muted text: Bright readable light silver */
      [data-bs-theme='dark'] .text-muted,
      [data-bs-theme='dark'] small.text-muted,
      [data-bs-theme='dark'] .text-secondary,
      [data-bs-theme='dark'] .stat-label {
        color: #CBD5E1 !important;
      }

      /* 3. High-Contrast Dark Mode Alert Banners & Callouts */
      [data-bs-theme='dark'] .alert {
        color: #F8FAFC !important;
        border-width: 1px !important;
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.35) !important;
      }
      [data-bs-theme='dark'] .alert p,
      [data-bs-theme='dark'] .alert div,
      [data-bs-theme='dark'] .alert span:not(.badge) {
        color: #F8FAFC !important;
      }
      [data-bs-theme='dark'] .alert b,
      [data-bs-theme='dark'] .alert strong,
      [data-bs-theme='dark'] .alert .alert-heading {
        color: #FFFFFF !important;
      }

      /* Alert Secondary (e.g. Reproducible Analysis Callout Banner) */
      [data-bs-theme='dark'] .alert-secondary {
        background-color: #1E293B !important;
        border-color: #475569 !important;
        color: #F8FAFC !important;
      }
      [data-bs-theme='dark'] .alert-secondary b,
      [data-bs-theme='dark'] .alert-secondary strong {
        color: #FFFFFF !important;
      }
      [data-bs-theme='dark'] .alert-secondary span:not(.badge) {
        color: #E2E8F0 !important;
      }

      /* Alert Info */
      [data-bs-theme='dark'] .alert-info {
        background-color: #082F49 !important;
        border-color: #0284C7 !important;
        color: #F0F9FF !important;
      }
      [data-bs-theme='dark'] .alert-info b,
      [data-bs-theme='dark'] .alert-info strong,
      [data-bs-theme='dark'] .alert-info .alert-heading {
        color: #FFFFFF !important;
      }
      [data-bs-theme='dark'] .alert-info a,
      [data-bs-theme='dark'] .alert-info actionLink,
      [data-bs-theme='dark'] .alert-info .action-button {
        color: #38BDF8 !important;
        font-weight: 700 !important;
      }

      /* Alert Primary & Indigo */
      [data-bs-theme='dark'] .alert-primary,
      [data-bs-theme='dark'] .alert-indigo {
        background-color: #1E1B4B !important;
        border-color: #6366F1 !important;
        color: #EEF2FF !important;
      }
      [data-bs-theme='dark'] .alert-primary b,
      [data-bs-theme='dark'] .alert-primary strong,
      [data-bs-theme='dark'] .alert-indigo b,
      [data-bs-theme='dark'] .alert-indigo strong {
        color: #FFFFFF !important;
      }
      [data-bs-theme='dark'] .alert-primary a,
      [data-bs-theme='dark'] .alert-primary actionLink,
      [data-bs-theme='dark'] .alert-primary .action-button {
        color: #A5B4FC !important;
        font-weight: 700 !important;
      }

      /* Alert Success */
      [data-bs-theme='dark'] .alert-success {
        background-color: #064E3B !important;
        border-color: #059669 !important;
        color: #ECFDF5 !important;
      }
      [data-bs-theme='dark'] .alert-success b,
      [data-bs-theme='dark'] .alert-success strong {
        color: #FFFFFF !important;
      }

      /* Alert Warning */
      [data-bs-theme='dark'] .alert-warning {
        background-color: #451A03 !important;
        border-color: #D97706 !important;
        color: #FEF3C7 !important;
      }
      [data-bs-theme='dark'] .alert-warning b,
      [data-bs-theme='dark'] .alert-warning strong {
        color: #FFFFFF !important;
      }

      /* Alert Danger */
      [data-bs-theme='dark'] .alert-danger {
        background-color: #450A0A !important;
        border-color: #DC2626 !important;
        color: #FEE2E2 !important;
      }
      [data-bs-theme='dark'] .alert-danger b,
      [data-bs-theme='dark'] .alert-danger strong {
        color: #FFFFFF !important;
      }

      /* 4. Attribute Selectors to Override Hardcoded Dark / Black Inline Text */
      [data-bs-theme='dark'] [style*='color: #0F172A'],
      [data-bs-theme='dark'] [style*='color: #0F172A'] *,
      [data-bs-theme='dark'] [style*='color:#0F172A'],
      [data-bs-theme='dark'] [style*='color:#0F172A'] *,
      [data-bs-theme='dark'] [style*='color: #0f172a'],
      [data-bs-theme='dark'] [style*='color: #0f172a'] *,
      [data-bs-theme='dark'] [style*='color:#0f172a'],
      [data-bs-theme='dark'] [style*='color:#0f172a'] *,
      [data-bs-theme='dark'] [style*='color: #1E293B'],
      [data-bs-theme='dark'] [style*='color: #1E293B'] *,
      [data-bs-theme='dark'] [style*='color:#1E293B'],
      [data-bs-theme='dark'] [style*='color:#1E293B'] *,
      [data-bs-theme='dark'] [style*='color: #1e293b'],
      [data-bs-theme='dark'] [style*='color: #1e293b'] *,
      [data-bs-theme='dark'] [style*='color:#1e293b'],
      [data-bs-theme='dark'] [style*='color:#1e293b'] *,
      [data-bs-theme='dark'] [style*='color: #334155'],
      [data-bs-theme='dark'] [style*='color: #334155'] *,
      [data-bs-theme='dark'] [style*='color:#334155'],
      [data-bs-theme='dark'] [style*='color:#334155'] *,
      [data-bs-theme='dark'] [style*='color: #475569'],
      [data-bs-theme='dark'] [style*='color: #475569'] *,
      [data-bs-theme='dark'] [style*='color:#475569'],
      [data-bs-theme='dark'] [style*='color:#475569'] *,
      [data-bs-theme='dark'] [style*='color: #111827'],
      [data-bs-theme='dark'] [style*='color: #111827'] *,
      [data-bs-theme='dark'] [style*='color:#111827'],
      [data-bs-theme='dark'] [style*='color:#111827'] *,
      [data-bs-theme='dark'] [style*='color: #000'],
      [data-bs-theme='dark'] [style*='color: #000'] *,
      [data-bs-theme='dark'] [style*='color:#000'],
      [data-bs-theme='dark'] [style*='color:#000'] *,
      [data-bs-theme='dark'] [style*='color: black'],
      [data-bs-theme='dark'] [style*='color: black'] * {
        color: #F8FAFC !important;
      }

      [data-bs-theme='dark'] [style*='color: #64748B'],
      [data-bs-theme='dark'] [style*='color: #64748B'] *,
      [data-bs-theme='dark'] [style*='color:#64748B'],
      [data-bs-theme='dark'] [style*='color:#64748B'] *,
      [data-bs-theme='dark'] [style*='color: #64748b'],
      [data-bs-theme='dark'] [style*='color: #64748b'] *,
      [data-bs-theme='dark'] [style*='color:#64748b'],
      [data-bs-theme='dark'] [style*='color:#64748b'] *,
      [data-bs-theme='dark'] [style*='color: #94A3B8'],
      [data-bs-theme='dark'] [style*='color: #94A3B8'] *,
      [data-bs-theme='dark'] [style*='color:#94A3B8'],
      [data-bs-theme='dark'] [style*='color:#94A3B8'] *,
      [data-bs-theme='dark'] [style*='color: #94a3b8'],
      [data-bs-theme='dark'] [style*='color: #94a3b8'] *,
      [data-bs-theme='dark'] [style*='color:#94a3b8'],
      [data-bs-theme='dark'] [style*='color:#94a3b8'] * {
        color: #CBD5E1 !important;
      }

      /* Hardcoded dark blues, navy, indigos become vibrant cyan / sky blue */
      [data-bs-theme='dark'] [style*='color: #1E3A8A'],
      [data-bs-theme='dark'] [style*='color: #1E3A8A'] *,
      [data-bs-theme='dark'] [style*='color:#1E3A8A'],
      [data-bs-theme='dark'] [style*='color:#1E3A8A'] *,
      [data-bs-theme='dark'] [style*='color: #1e3a8a'],
      [data-bs-theme='dark'] [style*='color: #1e3a8a'] *,
      [data-bs-theme='dark'] [style*='color:#1e3a8a'],
      [data-bs-theme='dark'] [style*='color:#1e3a8a'] *,
      [data-bs-theme='dark'] [style*='color: #1B4F72'],
      [data-bs-theme='dark'] [style*='color: #1B4F72'] *,
      [data-bs-theme='dark'] [style*='color:#1B4F72'],
      [data-bs-theme='dark'] [style*='color:#1B4F72'] *,
      [data-bs-theme='dark'] [style*='color: #1b4f72'],
      [data-bs-theme='dark'] [style*='color: #1b4f72'] *,
      [data-bs-theme='dark'] [style*='color:#1b4f72'],
      [data-bs-theme='dark'] [style*='color:#1b4f72'] *,
      [data-bs-theme='dark'] [style*='color: #1E3A5F'],
      [data-bs-theme='dark'] [style*='color: #1E3A5F'] *,
      [data-bs-theme='dark'] [style*='color:#1E3A5F'],
      [data-bs-theme='dark'] [style*='color:#1E3A5F'] *,
      [data-bs-theme='dark'] [style*='color: #1e3a5f'],
      [data-bs-theme='dark'] [style*='color: #1e3a5f'] *,
      [data-bs-theme='dark'] [style*='color:#1e3a5f'],
      [data-bs-theme='dark'] [style*='color:#1e3a5f'] *,
      [data-bs-theme='dark'] [style*='color: #172554'],
      [data-bs-theme='dark'] [style*='color: #172554'] *,
      [data-bs-theme='dark'] [style*='color:#172554'],
      [data-bs-theme='dark'] [style*='color:#172554'] *,
      [data-bs-theme='dark'] [style*='color: #1e40af'],
      [data-bs-theme='dark'] [style*='color: #1e40af'] *,
      [data-bs-theme='dark'] [style*='color:#1e40af'],
      [data-bs-theme='dark'] [style*='color:#1e40af'] *,
      [data-bs-theme='dark'] [style*='color: #1E40AF'],
      [data-bs-theme='dark'] [style*='color: #1E40AF'] *,
      [data-bs-theme='dark'] [style*='color:#1E40AF'],
      [data-bs-theme='dark'] [style*='color:#1E40AF'] *,
      [data-bs-theme='dark'] [style*='color: #0369A1'],
      [data-bs-theme='dark'] [style*='color: #0369A1'] *,
      [data-bs-theme='dark'] [style*='color:#0369A1'],
      [data-bs-theme='dark'] [style*='color:#0369A1'] *,
      [data-bs-theme='dark'] [style*='color: #0369a1'],
      [data-bs-theme='dark'] [style*='color: #0369a1'] *,
      [data-bs-theme='dark'] [style*='color:#0369a1'],
      [data-bs-theme='dark'] [style*='color:#0369a1'] *,
      [data-bs-theme='dark'] [style*='color: #0C4A6E'],
      [data-bs-theme='dark'] [style*='color: #0C4A6E'] *,
      [data-bs-theme='dark'] [style*='color:#0C4A6E'],
      [data-bs-theme='dark'] [style*='color:#0C4A6E'] *,
      [data-bs-theme='dark'] [style*='color: #0c4a6e'],
      [data-bs-theme='dark'] [style*='color: #0c4a6e'] *,
      [data-bs-theme='dark'] [style*='color:#0c4a6e'],
      [data-bs-theme='dark'] [style*='color:#0c4a6e'] *,
      [data-bs-theme='dark'] [style*='color: #082F49'],
      [data-bs-theme='dark'] [style*='color: #082F49'] *,
      [data-bs-theme='dark'] [style*='color:#082F49'],
      [data-bs-theme='dark'] [style*='color:#082F49'] *,
      [data-bs-theme='dark'] [style*='color: #2563EB'],
      [data-bs-theme='dark'] [style*='color: #2563EB'] *,
      [data-bs-theme='dark'] [style*='color:#2563EB'],
      [data-bs-theme='dark'] [style*='color:#2563EB'] *,
      [data-bs-theme='dark'] [style*='color: #3B82F6'],
      [data-bs-theme='dark'] [style*='color: #3B82F6'] *,
      [data-bs-theme='dark'] [style*='color:#3B82F6'],
      [data-bs-theme='dark'] [style*='color:#3B82F6'] * {
        color: #38BDF8 !important;
      }

      /* Hardcoded deep indigos become radiant soft indigo / violet */
      [data-bs-theme='dark'] [style*='color: #312E81'],
      [data-bs-theme='dark'] [style*='color: #312E81'] *,
      [data-bs-theme='dark'] [style*='color:#312E81'],
      [data-bs-theme='dark'] [style*='color:#312E81'] *,
      [data-bs-theme='dark'] [style*='color: #312e81'],
      [data-bs-theme='dark'] [style*='color: #312e81'] *,
      [data-bs-theme='dark'] [style*='color:#312e81'],
      [data-bs-theme='dark'] [style*='color:#312e81'] *,
      [data-bs-theme='dark'] [style*='color: #3730A3'],
      [data-bs-theme='dark'] [style*='color: #3730A3'] *,
      [data-bs-theme='dark'] [style*='color:#3730A3'],
      [data-bs-theme='dark'] [style*='color:#3730A3'] *,
      [data-bs-theme='dark'] [style*='color: #3730a3'],
      [data-bs-theme='dark'] [style*='color: #3730a3'] *,
      [data-bs-theme='dark'] [style*='color:#3730a3'],
      [data-bs-theme='dark'] [style*='color:#3730a3'] *,
      [data-bs-theme='dark'] [style*='color: #1E1B4B'],
      [data-bs-theme='dark'] [style*='color: #1E1B4B'] *,
      [data-bs-theme='dark'] [style*='color:#1E1B4B'],
      [data-bs-theme='dark'] [style*='color:#1E1B4B'] *,
      [data-bs-theme='dark'] [style*='color: #1e1b4b'],
      [data-bs-theme='dark'] [style*='color: #1e1b4b'] *,
      [data-bs-theme='dark'] [style*='color:#1e1b4b'],
      [data-bs-theme='dark'] [style*='color:#1e1b4b'] *,
      [data-bs-theme='dark'] [style*='color: #4338CA'],
      [data-bs-theme='dark'] [style*='color: #4338CA'] *,
      [data-bs-theme='dark'] [style*='color:#4338CA'],
      [data-bs-theme='dark'] [style*='color:#4338CA'] *,
      [data-bs-theme='dark'] [style*='color: #4338ca'],
      [data-bs-theme='dark'] [style*='color: #4338ca'] *,
      [data-bs-theme='dark'] [style*='color:#4338ca'],
      [data-bs-theme='dark'] [style*='color:#4338ca'] *,
      [data-bs-theme='dark'] [style*='color: #4F46E5'],
      [data-bs-theme='dark'] [style*='color: #4F46E5'] *,
      [data-bs-theme='dark'] [style*='color:#4F46E5'],
      [data-bs-theme='dark'] [style*='color:#4F46E5'] *,
      [data-bs-theme='dark'] [style*='color: #4f46e5'],
      [data-bs-theme='dark'] [style*='color: #4f46e5'] *,
      [data-bs-theme='dark'] [style*='color:#4f46e5'],
      [data-bs-theme='dark'] [style*='color:#4f46e5'] * {
        color: #A5B4FC !important;
      }

      /* Hardcoded dark greens become radiant emerald */
      [data-bs-theme='dark'] [style*='color: #166534'],
      [data-bs-theme='dark'] [style*='color: #166534'] *,
      [data-bs-theme='dark'] [style*='color:#166534'],
      [data-bs-theme='dark'] [style*='color:#166534'] *,
      [data-bs-theme='dark'] [style*='color: #14532D'],
      [data-bs-theme='dark'] [style*='color: #14532D'] *,
      [data-bs-theme='dark'] [style*='color:#14532D'],
      [data-bs-theme='dark'] [style*='color:#14532D'] *,
      [data-bs-theme='dark'] [style*='color: #14532d'],
      [data-bs-theme='dark'] [style*='color: #14532d'] *,
      [data-bs-theme='dark'] [style*='color:#14532d'],
      [data-bs-theme='dark'] [style*='color:#14532d'] *,
      [data-bs-theme='dark'] [style*='color: #15803D'],
      [data-bs-theme='dark'] [style*='color: #15803D'] *,
      [data-bs-theme='dark'] [style*='color:#15803D'],
      [data-bs-theme='dark'] [style*='color:#15803D'] *,
      [data-bs-theme='dark'] [style*='color: #15803d'],
      [data-bs-theme='dark'] [style*='color: #15803d'] *,
      [data-bs-theme='dark'] [style*='color:#15803d'],
      [data-bs-theme='dark'] [style*='color:#15803d'] *,
      [data-bs-theme='dark'] [style*='color: #16A34A'],
      [data-bs-theme='dark'] [style*='color: #16A34A'] *,
      [data-bs-theme='dark'] [style*='color:#16A34A'],
      [data-bs-theme='dark'] [style*='color:#16A34A'] *,
      [data-bs-theme='dark'] [style*='color: #16a34a'],
      [data-bs-theme='dark'] [style*='color: #16a34a'] *,
      [data-bs-theme='dark'] [style*='color:#16a34a'],
      [data-bs-theme='dark'] [style*='color:#16a34a'] *,
      [data-bs-theme='dark'] [style*='color: #0F766E'],
      [data-bs-theme='dark'] [style*='color: #0F766E'] *,
      [data-bs-theme='dark'] [style*='color:#0F766E'],
      [data-bs-theme='dark'] [style*='color:#0F766E'] *,
      [data-bs-theme='dark'] [style*='color: #0f766e'],
      [data-bs-theme='dark'] [style*='color: #0f766e'] *,
      [data-bs-theme='dark'] [style*='color:#0f766e'],
      [data-bs-theme='dark'] [style*='color:#0f766e'] *,
      [data-bs-theme='dark'] [style*='color: #0D9488'],
      [data-bs-theme='dark'] [style*='color: #0D9488'] *,
      [data-bs-theme='dark'] [style*='color:#0D9488'],
      [data-bs-theme='dark'] [style*='color:#0D9488'] *,
      [data-bs-theme='dark'] [style*='color: #0d9488'],
      [data-bs-theme='dark'] [style*='color: #0d9488'] *,
      [data-bs-theme='dark'] [style*='color:#0d9488'],
      [data-bs-theme='dark'] [style*='color:#0d9488'] * {
        color: #4ADE80 !important;
      }

      /* Hardcoded amber/browns become soft light gold */
      [data-bs-theme='dark'] [style*='color: #92400E'],
      [data-bs-theme='dark'] [style*='color: #92400E'] *,
      [data-bs-theme='dark'] [style*='color:#92400E'],
      [data-bs-theme='dark'] [style*='color:#92400E'] *,
      [data-bs-theme='dark'] [style*='color: #92400e'],
      [data-bs-theme='dark'] [style*='color: #92400e'] *,
      [data-bs-theme='dark'] [style*='color:#92400e'],
      [data-bs-theme='dark'] [style*='color:#92400e'] *,
      [data-bs-theme='dark'] [style*='color: #78350F'],
      [data-bs-theme='dark'] [style*='color: #78350F'] *,
      [data-bs-theme='dark'] [style*='color:#78350F'],
      [data-bs-theme='dark'] [style*='color:#78350F'] *,
      [data-bs-theme='dark'] [style*='color: #78350f'],
      [data-bs-theme='dark'] [style*='color: #78350f'] *,
      [data-bs-theme='dark'] [style*='color:#78350f'],
      [data-bs-theme='dark'] [style*='color:#78350f'] * {
        color: #FCD34D !important;
      }

      /* Hardcoded dark reds become soft coral / rose */
      [data-bs-theme='dark'] [style*='color: #991B1B'],
      [data-bs-theme='dark'] [style*='color: #991B1B'] *,
      [data-bs-theme='dark'] [style*='color:#991B1B'],
      [data-bs-theme='dark'] [style*='color:#991B1B'] *,
      [data-bs-theme='dark'] [style*='color: #991b1b'],
      [data-bs-theme='dark'] [style*='color: #991b1b'] *,
      [data-bs-theme='dark'] [style*='color:#991b1b'],
      [data-bs-theme='dark'] [style*='color:#991b1b'] *,
      [data-bs-theme='dark'] [style*='color: #7F1D1D'],
      [data-bs-theme='dark'] [style*='color: #7F1D1D'] *,
      [data-bs-theme='dark'] [style*='color:#7F1D1D'],
      [data-bs-theme='dark'] [style*='color:#7F1D1D'] *,
      [data-bs-theme='dark'] [style*='color: #7f1d1d'],
      [data-bs-theme='dark'] [style*='color: #7f1d1d'] *,
      [data-bs-theme='dark'] [style*='color:#7f1d1d'],
      [data-bs-theme='dark'] [style*='color:#7f1d1d'] *,
      [data-bs-theme='dark'] [style*='color: #DC3545'],
      [data-bs-theme='dark'] [style*='color: #DC3545'] *,
      [data-bs-theme='dark'] [style*='color:#DC3545'],
      [data-bs-theme='dark'] [style*='color:#DC3545'] *,
      [data-bs-theme='dark'] [style*='color: #dc3545'],
      [data-bs-theme='dark'] [style*='color: #dc3545'] *,
      [data-bs-theme='dark'] [style*='color:#dc3545'],
      [data-bs-theme='dark'] [style*='color:#dc3545'] * {
        color: #F87171 !important;
      }

      /* 5. Transform Hardcoded Light Background Containers into Dark Slate */
      [data-bs-theme='dark'] [style*='background: #FFFFFF'],
      [data-bs-theme='dark'] [style*='background:#FFFFFF'],
      [data-bs-theme='dark'] [style*='background-color: #FFFFFF'],
      [data-bs-theme='dark'] [style*='background-color:#FFFFFF'],
      [data-bs-theme='dark'] [style*='background: white'],
      [data-bs-theme='dark'] [style*='background-color: white'],
      [data-bs-theme='dark'] [style*='background: #F8FAFC'],
      [data-bs-theme='dark'] [style*='background:#F8FAFC'],
      [data-bs-theme='dark'] [style*='background-color: #F8FAFC'],
      [data-bs-theme='dark'] [style*='background-color:#F8FAFC'],
      [data-bs-theme='dark'] [style*='background: #EEF2F6'],
      [data-bs-theme='dark'] [style*='background:#EEF2F6'],
      [data-bs-theme='dark'] .bg-white,
      [data-bs-theme='dark'] .bg-light {
        background-color: #1E293B !important;
        color: #F8FAFC !important;
      }

      [data-bs-theme='dark'] [style*='background: #EFF6FF'],
      [data-bs-theme='dark'] [style*='background:#EFF6FF'],
      [data-bs-theme='dark'] [style*='background: #F0F9FF'],
      [data-bs-theme='dark'] [style*='background:#F0F9FF'] {
        background: #0C4A6E !important;
        color: #E0F2FE !important;
      }

      [data-bs-theme='dark'] [style*='background: #EEF2FF'],
      [data-bs-theme='dark'] [style*='background:#EEF2FF'],
      [data-bs-theme='dark'] [style*='background-color: #EEF2FF'],
      [data-bs-theme='dark'] [style*='background-color:#EEF2FF'] {
        background: #1E1B4B !important;
        color: #E0E7FF !important;
      }

      [data-bs-theme='dark'] [style*='background: #F0FDF4'],
      [data-bs-theme='dark'] [style*='background:#F0FDF4'] {
        background: #064E3B !important;
        color: #D1FAE5 !important;
      }

      [data-bs-theme='dark'] [style*='background: #FFFBEB'],
      [data-bs-theme='dark'] [style*='background:#FFFBEB'] {
        background: #451A03 !important;
        color: #FEF3C7 !important;
      }

      [data-bs-theme='dark'] [style*='background: #FEF2F2'],
      [data-bs-theme='dark'] [style*='background:#FEF2F2'] {
        background: #450A0A !important;
        color: #FEE2E2 !important;
      }

      /* Header Gradients in Dark Mode */
      [data-bs-theme='dark'] [style*='linear-gradient(135deg, #EFF6FF'] {
        background: linear-gradient(135deg, #0C4A6E 0%, #082F49 100%) !important;
        border-color: #0284C7 !important;
      }
      [data-bs-theme='dark'] [style*='linear-gradient(135deg, #EFF6FF'] h5,
      [data-bs-theme='dark'] [style*='linear-gradient(135deg, #EFF6FF'] p {
        color: #E0F2FE !important;
      }
      [data-bs-theme='dark'] [style*='linear-gradient(135deg, #EEF2FF'] {
        background: linear-gradient(135deg, #1E1B4B 0%, #311042 100%) !important;
        border-color: #4F46E5 !important;
      }
      [data-bs-theme='dark'] [style*='linear-gradient(135deg, #EEF2FF'] h5,
      [data-bs-theme='dark'] [style*='linear-gradient(135deg, #EEF2FF'] p {
        color: #E0E7FF !important;
      }

      /* Borders in Dark Mode */
      [data-bs-theme='dark'] .border,
      [data-bs-theme='dark'] [style*='border: 1px solid #E2E8F0'],
      [data-bs-theme='dark'] [style*='border:1px solid #E2E8F0'],
      [data-bs-theme='dark'] [style*='border: 1px solid #CBD5E1'],
      [data-bs-theme='dark'] [style*='border:1px solid #CBD5E1'],
      [data-bs-theme='dark'] [style*='border: 2px dashed #CBD5E1'],
      [data-bs-theme='dark'] [style*='border: 1px dashed #CBD5E1'] {
        border-color: #334155 !important;
      }

      /* 6. Form Controls, Inputs, File Uploads, Selects */
      [data-bs-theme='dark'] .form-control,
      [data-bs-theme='dark'] .form-select,
      [data-bs-theme='dark'] .selectize-input,
      [data-bs-theme='dark'] .selectize-dropdown {
        background-color: #1E293B !important;
        color: #F8FAFC !important;
        border-color: #475569 !important;
      }
      [data-bs-theme='dark'] .form-control:focus,
      [data-bs-theme='dark'] .form-select:focus,
      [data-bs-theme='dark'] .selectize-input.focus {
        background-color: #1E293B !important;
        color: #FFFFFF !important;
        border-color: #38BDF8 !important;
        box-shadow: 0 0 0 0.25rem rgba(56, 189, 248, 0.25) !important;
      }
      [data-bs-theme='dark'] .form-control::placeholder {
        color: #64748B !important;
      }
      [data-bs-theme='dark'] .form-check-input {
        background-color: #1E293B !important;
        border-color: #475569 !important;
      }
      [data-bs-theme='dark'] .form-check-input:checked {
        background-color: #0284C7 !important;
        border-color: #0284C7 !important;
      }
      [data-bs-theme='dark'] .form-check-label {
        color: #F8FAFC !important;
      }

      /* 7. Navigation Bar & Tabs */
      [data-bs-theme='dark'] .navbar {
        background-color: #0F172A !important;
        border-bottom: 1px solid #334155 !important;
      }
      [data-bs-theme='dark'] .navbar-brand {
        color: #FFFFFF !important;
      }
      [data-bs-theme='dark'] .navbar .nav-link {
        color: #CBD5E1 !important;
      }
      [data-bs-theme='dark'] .navbar .nav-link:hover,
      [data-bs-theme='dark'] .navbar .nav-link:focus {
        color: #38BDF8 !important;
      }
      [data-bs-theme='dark'] .navbar .nav-link.active {
        color: #38BDF8 !important;
        font-weight: 700 !important;
      }
      [data-bs-theme='dark'] .nav-pills .nav-link {
        color: #CBD5E1 !important;
        background-color: #1E293B !important;
        border: 1px solid #334155 !important;
      }
      [data-bs-theme='dark'] .nav-pills .nav-link.active {
        background-color: #0284C7 !important;
        color: #FFFFFF !important;
        border-color: #0284C7 !important;
      }

      /* 8. Cards, Placeholders, Containers */
      [data-bs-theme='dark'] .card {
        background-color: #0F172A !important;
        color: #F8FAFC !important;
        border-color: #334155 !important;
      }
      [data-bs-theme='dark'] .card-header {
        background-color: #1E293B !important;
        color: #FFFFFF !important;
        border-bottom-color: #334155 !important;
      }
      [data-bs-theme='dark'] .card-body {
        background-color: transparent !important;
        color: #F8FAFC !important;
      }
      [data-bs-theme='dark'] .stat-card {
        background: #1E293B !important;
        border-left-color: #38BDF8 !important;
        box-shadow: 0 2px 8px rgba(0,0,0,0.4) !important;
      }
      [data-bs-theme='dark'] .stat-number {
        color: #FFFFFF !important;
      }
      [data-bs-theme='dark'] .bubble-scroll-container {
        background: #0F172A !important;
        border-color: #334155 !important;
      }
      [data-bs-theme='dark'] .guide-step {
        background: #1E293B !important;
        border-color: #334155 !important;
        color: #F8FAFC !important;
      }
      [data-bs-theme='dark'] .info-card {
        background: #1E293B !important;
        border-color: #334155 !important;
        color: #F8FAFC !important;
      }
      [data-bs-theme='dark'] .info-card:hover {
        background: #24354D !important;
        border-color: #38BDF8 !important;
      }
      [data-bs-theme='dark'] .plot-placeholder-card {
        background: #1E293B !important;
        border: 2px dashed #475569 !important;
        box-shadow: 0 4px 14px rgba(0,0,0,0.3) !important;
      }
      [data-bs-theme='dark'] .plot-placeholder-card h5,
      [data-bs-theme='dark'] .plot-placeholder-card .placeholder-title {
        color: #FFFFFF !important;
      }
      [data-bs-theme='dark'] .plot-placeholder-card p,
      [data-bs-theme='dark'] .plot-placeholder-card .placeholder-desc {
        color: #CBD5E1 !important;
      }
      [data-bs-theme='dark'] .sim-input-card {
        background: #1E293B !important;
        border-color: #38BDF8 !important;
      }

      /* 9. Modals & Tables */
      [data-bs-theme='dark'] .modal-content {
        background-color: #0F172A !important;
        color: #F8FAFC !important;
        border-color: #334155 !important;
      }
      [data-bs-theme='dark'] .modal-header,
      [data-bs-theme='dark'] .modal-footer {
        border-color: #334155 !important;
        background-color: #1E293B !important;
      }
      [data-bs-theme='dark'] .table,
      [data-bs-theme='dark'] table.dataTable,
      [data-bs-theme='dark'] .table > :not(caption) > * > * {
        background-color: transparent !important;
        color: #F8FAFC !important;
        border-color: #334155 !important;
      }
      [data-bs-theme='dark'] .table-light,
      [data-bs-theme='dark'] thead.table-light th {
        background-color: #1E293B !important;
        color: #FFFFFF !important;
      }
      [data-bs-theme='dark'] pre,
      [data-bs-theme='dark'] code {
        background-color: #0F172A !important;
        color: #38BDF8 !important;
        border-color: #334155 !important;
      }

      /* DataTables Controls & Text in Dark Mode */
      [data-bs-theme='dark'] .dataTables_wrapper,
      [data-bs-theme='dark'] .dataTables_info,
      [data-bs-theme='dark'] .dataTables_length,
      [data-bs-theme='dark'] .dataTables_filter,
      [data-bs-theme='dark'] .dataTables_paginate,
      [data-bs-theme='dark'] .dataTables_filter label,
      [data-bs-theme='dark'] .dataTables_length label {
        color: #CBD5E1 !important;
      }
      [data-bs-theme='dark'] .dataTables_filter input,
      [data-bs-theme='dark'] .dataTables_length select {
        background-color: #1E293B !important;
        color: #F8FAFC !important;
        border: 1px solid #475569 !important;
      }
      [data-bs-theme='dark'] .dataTables_paginate .paginate_button {
        color: #CBD5E1 !important;
        background: #1E293B !important;
        border-color: #334155 !important;
      }
      [data-bs-theme='dark'] .dataTables_paginate .paginate_button.current,
      [data-bs-theme='dark'] .dataTables_paginate .paginate_button.current:hover {
        color: #FFFFFF !important;
        background: #0284C7 !important;
        border-color: #0284C7 !important;
      }
      [data-bs-theme='dark'] .dataTables_paginate .paginate_button:hover {
        color: #FFFFFF !important;
        background: #334155 !important;
        border-color: #475569 !important;
      }

      /* 10. Buttons & Interactive Links */
      [data-bs-theme='dark'] .btn-outline-primary {
        color: #38BDF8 !important;
        border-color: #38BDF8 !important;
        background-color: transparent !important;
      }
      [data-bs-theme='dark'] .btn-outline-primary:hover {
        background-color: #0284C7 !important;
        color: #FFFFFF !important;
      }
      [data-bs-theme='dark'] .btn-outline-secondary {
        color: #CBD5E1 !important;
        border-color: #64748B !important;
        background-color: transparent !important;
      }
      [data-bs-theme='dark'] .btn-outline-secondary:hover {
        background-color: #334155 !important;
        color: #FFFFFF !important;
      }
      [data-bs-theme='dark'] #btn_show_metadata_format_hint {
        background: #0F172A !important;
        color: #38BDF8 !important;
        border-color: #38BDF8 !important;
      }
      [data-bs-theme='dark'] #btn_show_metadata_format_multi_hint {
        background: #0F172A !important;
        color: #A5B4FC !important;
        border-color: #818CF8 !important;
      }
      [data-bs-theme='dark'] .badge.bg-white,
      [data-bs-theme='dark'] .badge.bg-light {
        background-color: #1E293B !important;
        color: #F8FAFC !important;
        border-color: #475569 !important;
      }
      [data-bs-theme='dark'] .pageview-box {
        background-color: rgba(30, 41, 59, 0.96) !important;
        border-color: #334155 !important;
        box-shadow: 0 4px 14px rgba(0,0,0,0.5) !important;
      }
      [data-bs-theme='dark'] .pageview-title {
        color: #94A3B8 !important;
      }
      [data-bs-theme='dark'] .pageview-count {
        color: #38BDF8 !important;
      }
      
      /* ===== LIVE LOGS TERMINAL CONSOLE STYLING ===== */
      .log-terminal-window {
        background: #090D16;
        border: 1px solid #1E293B;
        border-radius: 10px;
        box-shadow: 0 10px 30px rgba(0, 0, 0, 0.35);
        overflow: hidden;
      }
      .log-terminal-header {
        background: #111827;
        border-bottom: 1px solid #1E293B;
        padding: 9px 16px;
        display: flex;
        justify-content: space-between;
        align-items: center;
      }
      .terminal-dot {
        width: 11px;
        height: 11px;
        border-radius: 50%;
        display: inline-block;
      }
      .terminal-dot-red { background: #EF4444; }
      .terminal-dot-yellow { background: #F59E0B; }
      .terminal-dot-green { background: #10B981; }
      .log-terminal-body {
        background: #080C14;
        color: #E2E8F0;
        font-family: 'Consolas', 'Fira Code', 'Monaco', 'Courier New', monospace;
        font-size: 0.83rem;
        line-height: 1.6;
        padding: 16px 20px;
        max-height: 560px;
        min-height: 380px;
        overflow-y: auto;
        white-space: pre-wrap;
        word-break: break-word;
      }
      .log-entry-row {
        padding: 5px 0;
        border-bottom: 1px solid rgba(255, 255, 255, 0.04);
        display: flex;
        align-items: flex-start;
        gap: 10px;
      }
      .log-entry-row:hover {
        background: rgba(255, 255, 255, 0.025);
      }
      .log-time {
        color: #64748B;
        font-size: 0.77rem;
        min-width: 140px;
        flex-shrink: 0;
        user-select: none;
      }
      .log-lvl-badge {
        font-size: 0.70rem;
        font-weight: 700;
        padding: 2px 7px;
        border-radius: 4px;
        text-transform: uppercase;
        min-width: 72px;
        text-align: center;
        flex-shrink: 0;
        letter-spacing: 0.5px;
        user-select: none;
      }
      .log-lvl-ERROR {
        background: rgba(239, 68, 68, 0.22);
        color: #FCA5A5;
        border: 1px solid #DC2626;
      }
      .log-lvl-WARNING {
        background: rgba(245, 158, 11, 0.22);
        color: #FDE68A;
        border: 1px solid #D97706;
      }
      .log-lvl-NOTE {
        background: rgba(139, 92, 246, 0.22);
        color: #DDD6FE;
        border: 1px solid #7C3AED;
      }
      .log-lvl-INFO {
        background: rgba(56, 189, 248, 0.22);
        color: #7DD3FC;
        border: 1px solid #0284C7;
      }
      .log-lvl-SUCCESS {
        background: rgba(16, 185, 129, 0.22);
        color: #6EE7B7;
        border: 1px solid #059669;
      }
      .log-comp {
        color: #A5B4FC;
        font-weight: 600;
        min-width: 95px;
        flex-shrink: 0;
        font-size: 0.79rem;
      }
      .log-text {
        color: #F1F5F9;
        flex-grow: 1;
      }
      .log-details-block {
        margin-top: 5px;
        padding: 6px 10px;
        background: rgba(0, 0, 0, 0.4);
        border-left: 2px solid #EF4444;
        border-radius: 4px;
        color: #CBD5E1;
        font-size: 0.76rem;
      }
      @keyframes pulseLive {
        0% { opacity: 1; transform: scale(1); }
        50% { opacity: 0.35; transform: scale(0.85); }
        100% { opacity: 1; transform: scale(1); }
      }

      /* ===== PREVENT FILE UPLOAD AUTO-SCROLL TO TOP (SHINY INPUT FIX) ===== */
      html, body {
        scroll-behavior: auto !important;
      }
      .btn-file {
        position: relative !important;
        overflow: hidden !important;
        display: inline-block !important;
      }
      .btn-file input[type='file'],
      input[type='file'].shiny-input-file,
      .shiny-bound-input[type='file'] {
        position: absolute !important;
        top: 0 !important;
        left: 0 !important;
        right: 0 !important;
        bottom: 0 !important;
        width: 100% !important;
        height: 100% !important;
        min-width: 100% !important;
        min-height: 100% !important;
        margin: 0 !important;
        padding: 0 !important;
        opacity: 0 !important;
        cursor: pointer !important;
        display: block !important;
        visibility: visible !important;
        z-index: 9999 !important;
      }
      .input-group-btn,
      .input-group-prepend,
      .input-group-append {
        position: relative !important;
      }
    ")),
    tags$script(HTML("
      /* Proactive DOM Style Neuter & Continuous Scroll-Position Lock for File Inputs */
      (function() {
        // 1. Proactively neuter Shiny's inline style on file inputs
        function neuterFileInput(el) {
          if (!el) return;
          el.style.setProperty('position', 'absolute', 'important');
          el.style.setProperty('top', '0px', 'important');
          el.style.setProperty('left', '0px', 'important');
          el.style.setProperty('right', '0px', 'important');
          el.style.setProperty('bottom', '0px', 'important');
          el.style.setProperty('width', '100%', 'important');
          el.style.setProperty('height', '100%', 'important');
          el.style.setProperty('min-width', '100%', 'important');
          el.style.setProperty('min-height', '100%', 'important');
          el.style.setProperty('margin', '0px', 'important');
          el.style.setProperty('padding', '0px', 'important');
          el.style.setProperty('opacity', '0', 'important');
          el.style.setProperty('cursor', 'pointer', 'important');
          el.style.setProperty('display', 'block', 'important');
          el.style.setProperty('visibility', 'visible', 'important');
          el.style.setProperty('z-index', '9999', 'important');
          var p = el.parentElement;
          if (p) {
            p.style.setProperty('position', 'relative', 'important');
            p.style.setProperty('overflow', 'hidden', 'important');
          }
        }

        function neuterAllFileInputs() {
          var inputs = document.querySelectorAll('input[type=\"file\"]');
          for (var i = 0; i < inputs.length; i++) {
            neuterFileInput(inputs[i]);
          }
        }

        // Run immediately and observe DOM changes
        neuterAllFileInputs();
        if (window.MutationObserver) {
          var observer = new MutationObserver(function() {
            neuterAllFileInputs();
          });
          observer.observe(document.documentElement, {
            childList: true,
            subtree: true,
            attributes: true,
            attributeFilter: ['style', 'class']
          });
        }

        $(document).on('shiny:connected shiny:sessioninitialized shiny:value shiny:idle shiny:inputchanged', function() {
          neuterAllFileInputs();
        });

        // 2. High-precision scroll state tracking & continuous enforcement
        var activeScrollMap = [];
        var isLockActive = false;

        function recordCurrentScroll() {
          activeScrollMap = [];
          // Save window / document scroll
          var winY = window.pageYOffset || (document.documentElement ? document.documentElement.scrollTop : 0) || (document.body ? document.body.scrollTop : 0) || 0;
          var winX = window.pageXOffset || (document.documentElement ? document.documentElement.scrollLeft : 0) || (document.body ? document.body.scrollLeft : 0) || 0;
          activeScrollMap.push({ target: 'window', top: winY, left: winX });

          // Save all currently scrolled containers (.main, sidebar, layout-sidebar, etc.)
          var allNodes = document.querySelectorAll('*');
          for (var i = 0; i < allNodes.length; i++) {
            var n = allNodes[i];
            if (n.scrollTop > 0 || n.scrollLeft > 0) {
              activeScrollMap.push({ target: n, top: n.scrollTop, left: n.scrollLeft });
            }
          }
        }

        function applyRecordedScroll() {
          if (!activeScrollMap || activeScrollMap.length === 0) return;
          for (var i = 0; i < activeScrollMap.length; i++) {
            var item = activeScrollMap[i];
            if (item.target === 'window') {
              if (item.top > 0 || item.left > 0) {
                window.scrollTo({ top: item.top, left: item.left, behavior: 'instant' });
                if (document.documentElement) {
                  document.documentElement.scrollTop = item.top;
                  document.documentElement.scrollLeft = item.left;
                }
                if (document.body) {
                  document.body.scrollTop = item.top;
                  document.body.scrollLeft = item.left;
                }
              }
            } else if (item.target && item.target.isConnected) {
              item.target.scrollTop = item.top;
              item.target.scrollLeft = item.left;
            }
          }
        }

        function lockScrollDuringUpload() {
          neuterAllFileInputs();
          recordCurrentScroll();
          isLockActive = true;
          
          var startTime = performance.now();
          function frameLock(currentTime) {
            if (!isLockActive) return;
            applyRecordedScroll();
            if (currentTime - startTime < 2500) {
              requestAnimationFrame(frameLock);
            } else {
              isLockActive = false;
            }
          }
          requestAnimationFrame(frameLock);

          var intervals = [0, 20, 50, 100, 200, 400, 700, 1000, 1500, 2000, 2500];
          for (var k = 0; k < intervals.length; k++) {
            setTimeout(applyRecordedScroll, intervals[k]);
          }
        }

        // Intercept pointer/mouse events in CAPTURE phase before browser default actions
        var eventTypes = ['pointerdown', 'mousedown', 'click'];
        for (var eIdx = 0; eIdx < eventTypes.length; eIdx++) {
          document.addEventListener(eventTypes[eIdx], function(e) {
            var t = e.target;
            if (t && (t.matches('input[type=\"file\"]') || t.closest('.btn-file, .input-group-btn, .shiny-input-container, label[for*=\"file\"]'))) {
              lockScrollDuringUpload();
            }
          }, true);
        }

        // On hover or pointer enter, preemptively neuter file inputs so they are never at -99999px
        var hoverEvents = ['pointerover', 'mouseover', 'mouseenter'];
        for (var hIdx = 0; hIdx < hoverEvents.length; hIdx++) {
          document.addEventListener(hoverEvents[hIdx], function(e) {
            var t = e.target;
            if (t && (t.matches('input[type=\"file\"]') || t.closest('.btn-file, .input-group-btn, .shiny-input-container'))) {
              neuterAllFileInputs();
            }
          }, true);
        }

        // On focusin in capture phase, prevent scroll jumping
        document.addEventListener('focusin', function(e) {
          var t = e.target;
          if (t && (t.matches('input[type=\"file\"]') || t.closest('.btn-file, .input-group-btn, .shiny-input-container'))) {
            applyRecordedScroll();
          }
        }, true);

        // When window regains focus (file chooser dismissed/closed) or input changes
        window.addEventListener('focus', function() {
          applyRecordedScroll();
          setTimeout(applyRecordedScroll, 50);
          setTimeout(applyRecordedScroll, 150);
        });

        document.addEventListener('change', function(e) {
          var t = e.target;
          if (t && t.matches && t.matches('input[type=\"file\"]')) {
            applyRecordedScroll();
            setTimeout(applyRecordedScroll, 50);
            setTimeout(applyRecordedScroll, 150);
            setTimeout(applyRecordedScroll, 300);
            setTimeout(applyRecordedScroll, 600);
          }
        }, true);
      })();
    ")),
    div(
      id = "app_busy_pill",
      class = "app-global-loader",
      tags$span(class = "spinner-border spinner-border-sm text-light", role = "status", `aria-hidden` = "true"),
      tags$span("Processing & Rendering...")
    )
  ),
  
  # ============================================================================
  # TAB 1: HOME
  # ============================================================================
  nav_panel(
    "Home",
    div(
      style = "max-width: 1240px; margin: 0 auto; padding: 10px 0 30px 0;",
      
      # 1. Clean Scientific Hero Section
      div(
        style = "background: linear-gradient(135deg, #0F172A 0%, #1E293B 100%); color: #F8FAFC; border-radius: 12px; padding: 34px 38px; margin-bottom: 24px; box-shadow: 0 4px 16px rgba(15, 23, 42, 0.10); border: 1px solid #334155;",
        div(
          style = "max-width: 920px;",
          div(
            style = "display: inline-flex; align-items: center; gap: 8px; background: rgba(56, 189, 248, 0.12); border: 1px solid rgba(56, 189, 248, 0.28); border-radius: 20px; padding: 4px 14px; margin-bottom: 14px;",
            icon("dna", style = "color: #38BDF8; font-size: 0.82rem;"),
            span("Simulation Evaluation • scRNA-seq • scATAC-seq • Multiomics", style = "color: #38BDF8; font-size: 0.78rem; font-weight: 600; letter-spacing: 0.04em; text-transform: uppercase;")
          ),
          h2(
            "scSimEval Studio",
            style = "font-weight: 800; font-size: 2.1rem; letter-spacing: -0.03em; color: #FFFFFF; margin-bottom: 10px;"
          ),
          p(
            "This application is designed for simulation evaluation of single-cell RNA-seq (scRNA-seq), single-cell ATAC-seq (scATAC-seq), and single-cell multiomics (scRNA-seq + scATAC-seq) datasets by comparing them against real reference datasets across 62 evaluation measures.",
            style = "font-size: 1.02rem; color: #CBD5E1; line-height: 1.6; margin-bottom: 0;"
          )
        )
      ),
      
      # 2. Minimalist KPI Metrics Ribbon
      fluidRow(
        column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #2563EB; padding: 14px 18px; border-radius: 8px; background: #FFFFFF; margin-bottom: 22px;",
                      div(class = "stat-number", style = "color: #0F172A; font-size: 1.85rem; font-weight: 800; font-family: monospace;", "62"),
                      div(class = "stat-label", style = "color: #64748B; font-size: 0.74rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 4px;", "Evaluation Measures"))),
        column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #0D9488; padding: 14px 18px; border-radius: 8px; background: #FFFFFF; margin-bottom: 22px;",
                      div(class = "stat-number", style = "color: #0F172A; font-size: 1.85rem; font-weight: 800; font-family: monospace;", "8"),
                      div(class = "stat-label", style = "color: #64748B; font-size: 0.74rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 4px;", "Biological Dimensions"))),
        column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #7C3AED; padding: 14px 18px; border-radius: 8px; background: #FFFFFF; margin-bottom: 22px;",
                      div(class = "stat-number", style = "color: #0F172A; font-size: 1.85rem; font-weight: 800; font-family: monospace;", "3"),
                      div(class = "stat-label", style = "color: #64748B; font-size: 0.74rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 4px;", "Modalities (RNA / ATAC / Co-Assay)"))),
        column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #D97706; padding: 14px 18px; border-radius: 8px; background: #FFFFFF; margin-bottom: 22px;",
                      div(class = "stat-number", style = "color: #0F172A; font-size: 1.85rem; font-weight: 800; font-family: monospace;", "0.0 – 1.0"),
                      div(class = "stat-label", style = "color: #64748B; font-size: 0.74rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 4px;", "Normalized Fidelity Scale")))
      ),
      
      # 3. Two Balanced Scientific Cards
      fluidRow(
        column(
          7,
          card(
            style = "border: 1px solid #E2E8F0; border-radius: 10px; box-shadow: 0 1px 3px rgba(0,0,0,0.03); margin-bottom: 20px;",
            card_header(
              div(style = "display: flex; align-items: center; justify-content: space-between;",
                  span(icon("layer-group", class = "me-2 text-primary"), tags$b("Evaluation Framework (8 Biological Categories)")),
                  span(class = "badge bg-light text-secondary border", style = "font-size: 0.75rem;", "62 Measures Total"))
            ),
            card_body(
              style = "padding: 20px;",
              div(
                style = "display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin-bottom: 18px;",
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("I. Distributional Properties"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Gene/cell moments, library sparsity (14)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("II. Correlations & Zeros"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Gene-gene, cell-cell, kinetic noise (6)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("III. Cellular Structure"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Clustering ARI/NMI, silhouettes, k-NN (10)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("IV. Batch & Technical Mixing"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("kBET, CMS, LISI confounder alignment (7)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("V. Biological Downstream"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("DEG precision/recall, DV, pathway enrichment (15)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("VI. Lineage Dynamics"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Pseudotime Spearman/Kendall concordances (2)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("VII. Cross-Modal Coupling"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Peak-gene links, FOSCTTM, Match@1 (6)", style = "font-size: 0.75rem; color: #64748B;")),
                div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 9px 12px;",
                    div(tags$b("VIII. Scalability"), style = "font-size: 0.83rem; color: #1E293B; margin-bottom: 2px;"),
                    div("Compute duration & peak memory efficiency (2)", style = "font-size: 0.75rem; color: #64748B;"))
              ),
              p(
                "All metrics are standardized onto a direction-aware [0, 1] fidelity scale, where 1.0 represents perfect concordance with empirical biology. This enables unbiased ranking across simulator architectures across 62 evaluation measures.",
                style = "font-size: 0.85rem; color: #475569; margin-bottom: 0; line-height: 1.5;"
              )
            )
          )
        ),
        
        column(
          5,
          card(
            style = "border: 1px solid #E2E8F0; border-radius: 10px; box-shadow: 0 1px 3px rgba(0,0,0,0.03); margin-bottom: 20px;",
            card_header(
              span(icon("network-wired", class = "me-2 text-primary"), tags$b("Benchmarking Pipeline"))
            ),
            card_body(
              style = "padding: 20px;",
              div(
                style = "display: flex; flex-direction: column; gap: 14px;",
                div(
                  style = "display: flex; gap: 12px; align-items: flex-start;",
                  div(style = "flex-shrink: 0; width: 28px; height: 28px; border-radius: 50%; background: #EFF6FF; color: #2563EB; font-weight: 700; font-size: 0.82rem; display: flex; align-items: center; justify-content: center; border: 1px solid #BFDBFE;", "1"),
                  div(
                    div(tags$b("Data Ingestion"), style = "font-size: 0.86rem; color: #0F172A;"),
                    div("Load demo benchmark data (6 simulators) or upload count matrices (.rds, .csv, SingleCellExperiment, Seurat).", style = "font-size: 0.79rem; color: #64748B; line-height: 1.4;")
                  )
                ),
                div(
                  style = "display: flex; gap: 12px; align-items: flex-start;",
                  div(style = "flex-shrink: 0; width: 28px; height: 28px; border-radius: 50%; background: #EFF6FF; color: #2563EB; font-weight: 700; font-size: 0.82rem; display: flex; align-items: center; justify-content: center; border: 1px solid #BFDBFE;", "2"),
                  div(
                    div(tags$b("Evaluation & Ranking"), style = "font-size: 0.86rem; color: #0F172A;"),
                    div("Comprehensive evaluation of scRNA-seq, scATAC-seq, and multiomics simulations across 62 measures covering distributions, clustering, batch effects, and cross-modality coupling.", style = "font-size: 0.79rem; color: #64748B; line-height: 1.4;")
                  )
                ),
                div(
                  style = "display: flex; gap: 12px; align-items: flex-start;",
                  div(style = "flex-shrink: 0; width: 28px; height: 28px; border-radius: 50%; background: #EFF6FF; color: #2563EB; font-weight: 700; font-size: 0.82rem; display: flex; align-items: center; justify-content: center; border: 1px solid #BFDBFE;", "3"),
                  div(
                    div(tags$b("Visual Analytics & Export"), style = "font-size: 0.86rem; color: #0F172A;"),
                    div("Explore the comparative bubble matrix, diagnostic QC plots, ordinations, and export high-res figures and CSV/Excel reports.", style = "font-size: 0.79rem; color: #64748B; line-height: 1.4;")
                  )
                )
              ),
              hr(style = "margin: 18px 0 12px 0; border-color: #E2E8F0;"),
              div(
                style = "display: flex; justify-content: space-between; align-items: center;",
                span("Ready to evaluate?", style = "font-size: 0.82rem; color: #64748B; font-weight: 500;"),
                actionButton("btn_go_data_inline", "Go to Data Hub →", class = "btn btn-sm btn-outline-primary", style = "font-size: 0.8rem; font-weight: 600;")
              )
            )
          )
        )
      ),
      
      # 4. Minimalist Academic Footer Note
      div(
        style = "margin-top: 10px; padding: 14px 18px; background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 8px; display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 10px;",
        div(
          style = "font-size: 0.8rem; color: #64748B;",
          tags$b("scSimEval"), " • Single-Cell & Multiomics Simulation Benchmarking Studio • Open-Source R Package"
        ),
        div(
          style = "font-size: 0.8rem; color: #64748B;",
          "Local Studio: ", tags$code("scSimEval::launch_scSimEval_app()", style = "background: #E2E8F0; color: #334155; padding: 2px 6px; border-radius: 4px; font-size: 0.76rem;")
        )
      )
    )
  ),
  
  # ============================================================================
  # ============================================================================
  # TAB 2: DATA HUB (INGESTION & BENCHMARK REGISTRY)
  # ============================================================================
  nav_panel(
    "Data Hub",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        title = div(icon("database", class = "me-2"), "Input Controls"),
        
        selectInput(
          "opt_data_mode", tags$b("Choose Data Mode:"),
          choices = c(
            "1. scRNA-seq datasets" = "scrna",
            "2. scATAC-seq datasets" = "scatac",
            "3. Single-Cell Multiomics datasets" = "multiomics",
            "Upload Saved Results (.rds / .csv)" = "upload_bench"
          ),
          selected = "multiomics",
          width = "100%"
        ),
        hr(style = "margin: 10px 0;"),
        
        # Mode 1: scRNA-seq datasets
        conditionalPanel(
          condition = "input.opt_data_mode == 'scrna'",
          div(
            style = "background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 6px; padding: 10px; margin-bottom: 12px;",
            div(style = "font-weight: 700; color: #1E3A8A; font-size: 0.86rem; margin-bottom: 6px;", icon("flask", class = "me-1"), "Load Example scRNA-seq Data:"),
            radioButtons(
              "opt_example_scrna_dataset", NULL,
              choices = c(
                "Rice Bud (Paired, 6 Simulators)" = "rice_bud",
                "Soybean Root (Unpaired, 3 Simulators)" = "soybean_root"
              ),
              selected = "rice_bud"
            ),
            actionButton("btn_load_example_scrna", "Load Example scRNA-seq Data", class = "btn btn-primary btn-sm w-100", icon = icon("play"))
          ),
          div(
            style = "font-size: 0.78rem; color: #1E3A8A; background: #EFF6FF; border-left: 3px solid #2563EB; border-radius: 4px; padding: 8px 10px; margin-bottom: 12px;",
            div(style = "font-weight: 700; color: #1E3A8A; margin-bottom: 2px;", icon("layer-group", class = "me-1"), "Custom scRNA-seq Evaluation:"),
            div("Configure reference matrix, simulations, runtime, and RAM in the Main Panel.")
          ),
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 6px; margin-bottom: 4px;", "1. Benchmarking Scope"),
          numericInput("num_uni_datasets", "Number of Datasets:", value = 1, min = 1, max = 50, step = 1),
          numericInput("num_uni_sims", "Number of Simulators:", value = 1, min = 1, max = 10, step = 1),
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 10px; margin-bottom: 4px;", "2. Compute & Acceleration"),
          numericInput("num_uni_cores", "CPU Cores / Threads:", value = default_n_cores, min = 1, max = max_sys_cores, step = 1),
          span(class = "text-muted", style = "font-size: 0.74rem; display: block; margin-top: -6px; margin-bottom: 8px;", sprintf("System detected: %d CPU cores", max_sys_cores)),
          hr(style = "margin: 12px 0;"),
          checkboxInput("chk_append_uni", "Add to current benchmark results", value = FALSE),
          actionButton("btn_run_uni_eval", "Calculate scRNA-seq Metrics", class = "btn btn-primary w-100 mt-2", icon = icon("play")),
          actionButton("btn_inspect_uni_data_side", "Inspect Dataset Properties", class = "btn btn-outline-secondary w-100 mt-2 btn-sm", icon = icon("table-cells"))
        ),
        
        # Mode 2: scATAC-seq datasets
        conditionalPanel(
          condition = "input.opt_data_mode == 'scatac'",
          div(
            style = "background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 6px; padding: 10px; margin-bottom: 12px;",
            div(style = "font-weight: 700; color: #0284C7; font-size: 0.86rem; margin-bottom: 6px;", icon("flask", class = "me-1"), "Load Example scATAC-seq Data:"),
            radioButtons(
              "opt_example_scatac_dataset", NULL,
              choices = c(
                "Rice Bud (Paired, 6 Simulators)" = "rice_bud",
                "Soybean Root (Unpaired, 3 Simulators)" = "soybean_root"
              ),
              selected = "rice_bud"
            ),
            actionButton("btn_load_example_scatac", "Load Example scATAC-seq Data", class = "btn btn-info text-white btn-sm w-100", icon = icon("play"))
          ),
          div(
            style = "font-size: 0.78rem; color: #0C4A6E; background: #F0F9FF; border-left: 3px solid #0284C7; border-radius: 4px; padding: 8px 10px; margin-bottom: 12px;",
            div(style = "font-weight: 700; color: #0C4A6E; margin-bottom: 2px;", icon("chart-simple", class = "me-1"), "Custom scATAC-seq Evaluation:"),
            div("Configure reference peak accessibility matrix, simulations, runtime, and RAM in the Main Panel.")
          ),
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 6px; margin-bottom: 4px;", "1. Benchmarking Scope"),
          numericInput("num_uni_datasets_atac", "Number of Datasets:", value = 1, min = 1, max = 50, step = 1),
          numericInput("num_uni_sims_atac", "Number of Simulators:", value = 1, min = 1, max = 10, step = 1),
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 10px; margin-bottom: 4px;", "2. Compute & Acceleration"),
          numericInput("num_uni_cores_atac", "CPU Cores / Threads:", value = default_n_cores, min = 1, max = max_sys_cores, step = 1),
          span(class = "text-muted", style = "font-size: 0.74rem; display: block; margin-top: -6px; margin-bottom: 8px;", sprintf("System detected: %d CPU cores", max_sys_cores)),
          hr(style = "margin: 12px 0;"),
          checkboxInput("chk_append_uni_atac", "Add to current benchmark results", value = FALSE),
          actionButton("btn_run_uni_eval_atac", "Calculate scATAC-seq Metrics", class = "btn btn-info text-white w-100 mt-2", icon = icon("play")),
          actionButton("btn_inspect_uni_data_side_atac", "Inspect Dataset Properties", class = "btn btn-outline-secondary w-100 mt-2 btn-sm", icon = icon("table-cells"))
        ),
        
        # Mode 3: Multiomics (scRNA-seq + scATAC-seq)
        conditionalPanel(
          condition = "input.opt_data_mode == 'multiomics'",
          radioButtons(
            "opt_multi_pairing", tags$b("Assay Type:"),
            choices = c(
              "Paired (Same Cells, 10x / SHARE-seq)" = "paired",
              "Unpaired (Independent Cells)" = "unpaired"
            ),
            selected = "paired"
          ),
          uiOutput("ui_pairing_info_banner"),
          div(
            style = "background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 6px; padding: 10px; margin-bottom: 12px;",
            div(style = "font-weight: 700; color: #312E81; font-size: 0.86rem; margin-bottom: 6px;", icon("flask", class = "me-1"), "Load Example Multiomics Data:"),
            conditionalPanel(
              condition = "input.opt_multi_pairing == 'paired'",
              p("Rice Bud (Oryza sativa) paired multiomics dataset (100 cells, 200 genes, 200 peaks, 8 cell types, 2 batches) evaluated across 6 simulators (scMultiSim, scDesign3, MOSim, scCross, scMoMtF, Matilda).",
                style = "font-size: 0.78rem; color: #475569; margin-bottom: 8px; line-height: 1.4;"),
              actionButton("btn_load_example_paired_multi", "Load Rice Bud Example (Paired, 6 Simulators)", class = "btn btn-primary btn-sm w-100", icon = icon("play"))
            ),
            conditionalPanel(
              condition = "input.opt_multi_pairing == 'unpaired'",
              p("Soybean Root (Glycine max) unpaired multiomics dataset (100 RNA cells, 100 ATAC cells, 200 genes, 200 peaks) evaluated across 3 simulators (scMultiSim, scDesign3, scCross).",
                style = "font-size: 0.78rem; color: #475569; margin-bottom: 8px; line-height: 1.4;"),
              actionButton("btn_load_example_unpaired_multi", "Load Soybean Root Example (Unpaired, 3 Simulators)", class = "btn btn-primary btn-sm w-100", icon = icon("play"))
            )
          ),
          div(
            style = "font-size: 0.78rem; color: #312E81; background: #EEF2FF; border-left: 3px solid #6366F1; border-radius: 4px; padding: 8px 10px; margin-bottom: 12px;",
            div(style = "font-weight: 700; color: #312E81; margin-bottom: 2px;", icon("dna", class = "me-1"), "Custom Multiomics Evaluation:"),
            div("Configure paired/unpaired RNA & ATAC matrices, simulations, runtime, and RAM in the Main Panel.")
          ),
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 6px; margin-bottom: 4px;", "1. Benchmarking Scope"),
          numericInput("num_multi_datasets", "Number of Datasets:", value = 1, min = 1, max = 30, step = 1),
          numericInput("num_multi_sims", "Number of Simulators:", value = 1, min = 1, max = 8, step = 1),
          tags$div(style = "font-size: 0.76rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 700; color: #475569; margin-top: 10px; margin-bottom: 4px;", "2. Compute & Acceleration"),
          numericInput("num_multi_cores", "CPU Cores / Threads:", value = default_n_cores, min = 1, max = max_sys_cores, step = 1),
          span(class = "text-muted", style = "font-size: 0.74rem; display: block; margin-top: -6px; margin-bottom: 8px;", sprintf("System detected: %d CPU cores", max_sys_cores)),
          hr(style = "margin: 12px 0;"),
          checkboxInput("chk_append_multi", "Add to current benchmark results", value = FALSE),
          actionButton("btn_run_multi_eval", "Calculate Multiomics Metrics", class = "btn btn-primary w-100 mt-2", icon = icon("dna")),
          actionButton("btn_inspect_multi_data_side", "Inspect Dataset Properties", class = "btn btn-outline-secondary w-100 mt-2 btn-sm", icon = icon("table-cells"))
        ),
        
        # Mode 4: Upload Saved Results (.rds / .csv)
        conditionalPanel(
          condition = "input.opt_data_mode == 'upload_bench'",
          p("Load a previously saved benchmark result file (.rds or .csv) from scSimEval.", style = "font-size: 0.82rem; color: #475569;"),
          fileInput("file_bench_upload", "Select Saved File(s) (.rds / .csv):", multiple = TRUE, accept = c(".rds", ".csv")),
          actionButton("btn_load_uploaded_bench", "Load Saved Results", class = "btn btn-secondary w-100", icon = icon("folder-open"))
        )
      ),
      uiOutput("ui_datahub_main_panel")
    )
  ),
  # ============================================================================
  # TAB 3: COMPARATIVE BUBBLE MATRIX (FLAGSHIP FIGURE)
  # ============================================================================
  nav_panel(
    "Comparative Bubble Matrix",
    layout_sidebar(
      sidebar = sidebar(
        width = 330,
        title = div(icon("sliders", class = "me-2 text-primary"), tags$strong("Display & Filters")),
        
        # Confirm & Plot Action Button
        div(
          class = "mb-3",
          actionButton(
            "btn_plot_bubble", "Confirm & Generate Plot",
            icon = icon("circle-dot"),
            class = "btn btn-primary w-100",
            style = "font-weight: 700; padding: 10px; font-size: 0.90rem; box-shadow: 0 2px 6px rgba(30,58,138,0.25);"
          )
        ),
        
        # Dataset Selection (All Datasets Average vs Individual Datasets)
        uiOutput("ui_bubble_dataset_selector"),
        uiOutput("ui_bubble_modality_selector"),
        hr(style = "margin: 10px 0;"),

        # 1. Filter Mode
        radioButtons(
          "rad_bubble_filter_mode", tags$b("Filter Mode:"),
          choices = c(
            "All 62 Measures"  = "all",
            "Custom Selection" = "custom"
          ),
          selected = "all"
        ),
        
        # 2. Custom Selection Panel
        conditionalPanel(
          condition = "input.rad_bubble_filter_mode == 'custom'",
          div(
            style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 8px; padding: 12px; margin-top: 4px; margin-bottom: 12px;",
            uiOutput("ui_custom_category_picker"),
            div(
              class = "d-flex justify-content-between align-items-center mt-2 mb-1",
              tags$b("2. Select Metrics:", style = "font-size: 0.88rem; color: #1E293B;"),
              div(
                actionLink("btn_custom_metrics_all", "All", style = "font-size: 0.78rem; margin-right: 6px; font-weight: 600; text-decoration: none;"),
                actionLink("btn_custom_metrics_none", "Clear", style = "font-size: 0.78rem; color: #dc3545; text-decoration: none;")
              )
            ),
            div(
              style = "max-height: 190px; overflow-y: auto; background: #FFFFFF; border: 1px solid #CBD5E1; border-radius: 6px; padding: 6px 10px;",
              uiOutput("ui_custom_metrics_checklist")
            ),
            div(
              style = "margin-top: 5px; font-size: 0.76rem; color: #64748B; text-align: right;",
              textOutput("txt_custom_metrics_count", inline = TRUE)
            )
          )
        ),
        
        hr(style = "margin: 10px 0;"),
        
        # 3. Simulator Methods Selection
        uiOutput("ui_bubble_method_picker"),
        
        hr(style = "margin: 10px 0;"),
        
        # 4. Collapsible Publication & Sizing Settings
        accordion(
          open = FALSE,
          accordion_panel(
            title = span(icon("sliders"), " Publication Title & Sizing"),
            value = "panel_dimensions",
            checkboxInput("chk_bubble_show_title", "Include Title & Subtitle in Plot", value = TRUE),
            textInput("txt_bubble_custom_title", "Custom Title (Optional):", value = "", placeholder = "e.g., Single-Cell Simulator Evaluation"),
            sliderInput("sld_bubble_width", "Matrix Width (px):", min = 800, max = 3200, value = 2200, step = 50),
            sliderInput("sld_bubble_height", "Matrix Height (px):", min = 350, max = 1200, value = 680, step = 20),
            p("Adjust plot dimensions or customize title for scientific publication.", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 0;")
          )
        ),
        
        hr(style = "margin: 12px 0;"),
        
        # 5. Clean download section
        div(
          class = "d-grid gap-2",
          downloadButton("download_bubble_jpeg", "Download JPEG (600 DPI)", class = "btn btn-primary btn-sm"),
          downloadButton("download_bubble_pdf", "Download Vector PDF", class = "btn btn-outline-secondary btn-sm"),
          downloadButton("download_bubble_r_script", "Download R Script", class = "btn btn-outline-primary btn-sm", icon = icon("file-code"))
        )
      ),
      
      navset_pill(
        id = "bubble_matrix_subtabs",
        selected = "Comparative Bubble Matrix Figure",
        
        # Subtab 1: Figure
        nav_panel(
          "Comparative Bubble Matrix Figure",
          card(
            card_header(
              div(
                class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
                span(icon("circle-dot", class = "me-2 text-primary"), tags$b("Simulator Comparison Bubble Matrix")),
                div(
                  class = "d-flex align-items-center gap-2",
                  uiOutput("ui_bubble_active_badge"),
                  downloadButton("download_bubble_r_script_header", "Download R Script", class = "btn btn-outline-primary btn-sm", icon = icon("file-code"))
                )
              )
            ),
            card_body(
              # Dataset Navigation Tabs (All Datasets Average vs Individual Datasets)
              uiOutput("ui_bubble_dataset_pills"),

              p("Each column represents a simulation method; each row represents an evaluation metric. Larger bubbles show higher scores (closer match to real data). Colors indicate metric categories.", style = "font-size: 0.84rem; color: #64748B; margin-bottom: 12px;"),
              div(
                class = "bubble-scroll-container",
                uiOutput("ui_bubble_plot_render")
              )
            )
          )
        ),
        
        # Subtab 2: Fidelity Leaderboard
        nav_panel(
          "Fidelity Leaderboard",
          card(
            card_header(
              class = "d-flex justify-content-between align-items-center py-2",
              tags$span(icon("ranking-star", class = "me-2 text-primary"), tags$b("Simulator Rankings")),
              downloadButton("download_leaderboard_csv", "Download Leaderboard (CSV)", class = "btn btn-sm btn-outline-primary")
            ),
            card_body(
              p(
                "Overall simulator rankings averaged across all evaluation metrics. Metrics where lower values are better (such as error, runtime, and memory) are inverted, and all scores are scaled between 0.00 and 1.00:",
                style = "font-size: 0.84rem; color: #64748B; margin-bottom: 12px;"
              ),
              DTOutput("table_leaderboard_dt")
            )
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 4: DIAGNOSTIC VISUALIZATIONS (8 SUB-PANELS)
  # ============================================================================
  nav_panel(
    "Visualizations",
    
    # Dataset selectors embedded directly within each visualization tab
    
    navset_pill(
      id = "viz_subtabs",
      
      # Sub-panel 1: Evaluation Summary
      nav_panel(
        "1. Evaluation Summary",
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              span(icon("chart-pie", class = "me-2 text-primary"), tags$b("Category Evaluation Summary")),
              div(
                class = "d-flex align-items-center gap-2",
                uiOutput("ui_viz_active_badge_sum"),
                downloadButton("download_sum_r_script_header", "Download R Script", class = "btn btn-outline-primary btn-sm", icon = icon("file-code")),
                span(class = "badge bg-light text-secondary border", "Radar & Bar Chart")
              )
            )
          ),
          card_body(
            fluidRow(
              column(3, uiOutput("ui_sum_dataset_selector")),
              column(2,
                     div(style = "margin-top: 24px;",
                         actionButton("btn_plot_sum", "Confirm & Plot", class = "btn btn-primary btn-sm w-100", icon = icon("chart-pie"), style = "font-weight: 600;")
                     )
              ),
              column(2, div(style = "margin-top: 28px;", checkboxInput("chk_sum_labels", "Show Score Labels", value = TRUE))),
              column(2, div(style = "margin-top: 28px;", checkboxInput("chk_sum_norm", "Normalize [0, 1]", value = TRUE))),
              column(3,
                     div(style = "margin-top: 20px;",
                         downloadButton("download_sum_jpeg", "JPEG (600 DPI)", class = "btn btn-sm btn-primary me-1 mb-1"),
                         downloadButton("download_sum_pdf", "PDF", class = "btn btn-sm btn-outline-secondary me-1 mb-1"),
                         downloadButton("download_sum_r_script", "R Script", class = "btn btn-sm btn-outline-primary mb-1", icon = icon("file-code"))
                     )
              )
            ),
            hr(),
            uiOutput("ui_plot_eval_summary")
          )
        )
      ),
      
      # Sub-panel 2: Distribution QC
      nav_panel(
        "2. Distribution QC",
        card(
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              span(icon("chart-line", class = "me-2 text-primary"), tags$b("Real vs Simulated Data Distributions")),
              div(
                class = "d-flex align-items-center gap-2",
                uiOutput("ui_viz_active_badge_dist"),
                downloadButton("download_dist_r_script_header", "Download R Script", class = "btn btn-outline-primary btn-sm", icon = icon("file-code")),
                span(class = "badge bg-light text-secondary border", "Density Curves")
              )
            )
          ),
          card_body(
            fluidRow(
              column(2,
                     actionButton("btn_plot_dist", "Confirm & Plot QC", class = "btn btn-primary btn-sm w-100", icon = icon("play"), style = "font-weight: 600;")
              ),
              column(2, uiOutput("ui_dist_dataset_selector")),
              column(2, uiOutput("ui_dist_modality_selector")),
              column(2, uiOutput("ui_dist_sim_selector")),
              column(2, selectInput("sel_dist_layout", "QC Layout:", choices = c("Comprehensive (14 Panels)" = "full", "Density Curves Only" = "density_only"), selected = "full")),
              column(2,
                     downloadButton("download_dist_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-1 mb-1"),
                     downloadButton("download_dist_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary me-1 mb-1"),
                     downloadButton("download_dist_r_script", "Download R Script", class = "btn btn-sm btn-outline-primary mb-1", icon = icon("file-code"))
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
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              span(icon("microchip", class = "me-2 text-primary"), tags$b("Runtime and Memory Usage")),
              div(
                class = "d-flex align-items-center gap-2",
                uiOutput("ui_viz_active_badge_scale"),
                downloadButton("download_scale_r_script_header", "Download R Script", class = "btn btn-outline-primary btn-sm", icon = icon("file-code")),
                span(class = "badge bg-light text-secondary border", "Computational Cost")
              )
            )
          ),
          card_body(
            fluidRow(
              column(3, uiOutput("ui_scale_dataset_selector")),
              column(3, uiOutput("ui_scale_modality_selector")),
              column(3,
                     selectInput(
                       "sel_scale_type", "Scalability View:",
                       choices = c(
                         "All 4 Panels Combined" = "composite",
                         "Runtime (Seconds)" = "runtime",
                         "Peak RAM Memory (MiB)" = "memory",
                         "Runtime vs Memory Trade-Off" = "tradeoff",
                         "Overall Resource Cost" = "cost"
                       ),
                       selected = "composite"
                     )
              ),
              column(3,
                     div(style = "margin-top: 24px;",
                         actionButton("btn_plot_scale", "Confirm & Plot", class = "btn btn-primary btn-sm w-100", icon = icon("microchip"), style = "font-weight: 600;")
                     )
              )
            ),
            fluidRow(
              column(12,
                     downloadButton("download_scale_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2 mb-1"),
                     downloadButton("download_scale_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary me-2 mb-1"),
                     downloadButton("download_scale_r_script", "Download R Script", class = "btn btn-sm btn-outline-primary mb-1", icon = icon("file-code"))
              )
            ),
            hr(),
            uiOutput("ui_plot_scale_bench")
          )
        )
      ),
      
      # Sub-panel 4: Metric Plots
      nav_panel(
        "4. Metric Plots",
        card(
          fill = FALSE,
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              span(icon("chart-column", class = "me-2 text-primary"), tags$b("Metric Plots by Category")),
              div(
                class = "d-flex align-items-center gap-2",
                uiOutput("ui_viz_active_badge_box"),
                downloadButton("download_box_r_script_header", "Download R Script", class = "btn btn-outline-primary btn-sm", icon = icon("file-code")),
                span(class = "badge bg-light text-secondary border", "Boxplots & Bars")
              )
            )
          ),
          card_body(
            fluidRow(
              column(8,
                     radioButtons(
                       "opt_box_view_mode", "View Mode:",
                       choices = c("View Individual Metric" = "individual", "View by Category Group" = "category"),
                       selected = "individual", inline = TRUE
                     )
              ),
              column(4,
                     actionButton("btn_plot_box", "Confirm & Plot Metric", class = "btn btn-primary btn-sm w-100", icon = icon("chart-column"), style = "font-weight: 600;")
              )
            ),
            fluidRow(
              conditionalPanel(
                condition = "input.opt_box_view_mode == 'individual'",
                column(2, uiOutput("ui_box_dataset_selector")),
                column(2, uiOutput("ui_box_modality_selector")),
                column(3,
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
                column(3, uiOutput("ui_box_metric_picker")),
                column(2,
                       selectInput(
                         "sel_box_indiv_score", "3. Score Type:",
                         choices = c("Normalized Score [0, 1]" = "normalized", "Raw Metric Value" = "raw"),
                         selected = "normalized"
                       )
                )
              ),
              conditionalPanel(
                condition = "input.opt_box_view_mode == 'category'",
                column(3, uiOutput("ui_box_dataset_selector_cat")),
                column(3, uiOutput("ui_box_modality_selector_cat")),
                column(3,
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
                column(3, selectInput("sel_box_score_type", "Score Type:", choices = c("Normalized Score [0, 1]" = "normalized", "Raw Value" = "raw"), selected = "normalized"))
              )
            ),
            fluidRow(
              column(12, style = "margin-top: 10px;",
                     downloadButton("download_box_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2 mb-1"),
                     downloadButton("download_box_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary me-2 mb-1"),
                     downloadButton("download_box_r_script", "Download R Script", class = "btn btn-sm btn-outline-primary mb-1", icon = icon("file-code"))
              )
            ),
            hr(),
            uiOutput("ui_plot_metric_boxes")
          )
        )
      ),
      
      # Sub-panel 5: Metric Heatmap
      nav_panel(
        "5. Metric Heatmap",
        card(
          fill = FALSE,
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              span(icon("table-cells", class = "me-2 text-primary"), tags$b("Performance Heatmap")),
              div(
                class = "d-flex align-items-center gap-2",
                uiOutput("ui_viz_active_badge_heat"),
                downloadButton("download_heat_r_script_header", "Download R Script", class = "btn btn-outline-primary btn-sm", icon = icon("file-code")),
                span(class = "badge bg-light text-secondary border", "Overview Table")
              )
            )
          ),
          card_body(
            fluidRow(
              column(3, uiOutput("ui_heat_dataset_selector")),
              column(2, uiOutput("ui_heat_modality_selector")),
              column(3,
                     div(style = "margin-bottom: 8px;",
                       actionButton("btn_plot_heat", "Confirm & Plot Heatmap", class = "btn btn-primary btn-sm w-100", icon = icon("table-cells"), style = "font-weight: 600;")
                     ),
                     selectInput(
                       "sel_heat_cat", "Category Filter:",
                       choices = c(
                         "All 8 Categories Combined" = "all",
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
                     selectInput(
                       "sel_heat_score_type", "Score Scale:",
                       choices = c("Normalized Fidelity [0, 1]" = "relative", "Raw Metric Values" = "raw"),
                       selected = "relative"
                     )
              ),
              column(2,
                     sliderInput("sld_heat_height", "Heatmap Height (px):", min = 400, max = 1500, value = 1100, step = 50)
              ),
              column(2,
                     sliderInput("sld_heat_width", "Heatmap Width (px):", min = 600, max = 1300, value = 880, step = 20)
              )
            ),
            fluidRow(
              column(12, style = "margin-top: 10px;",
                     downloadButton("download_heat_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2 mb-1"),
                     downloadButton("download_heat_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary me-2 mb-1"),
                     downloadButton("download_heat_r_script", "Download R Script", class = "btn btn-sm btn-outline-primary mb-1", icon = icon("file-code"))
              )
            ),
            hr(),
            div(
              class = "bubble-scroll-container",
              style = "text-align: center; overflow-x: auto; padding: 10px;",
              uiOutput("ui_plot_metric_heat")
            )
          )
        )
      ),
      
      # Sub-panel 6: PCA Ordination
      nav_panel(
        "6. PCA Ordination",
        card(
          fill = FALSE,
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              span(icon("compass", class = "me-2 text-primary"), tags$b("PCA Plot of Simulators")),
              div(
                class = "d-flex align-items-center gap-2",
                uiOutput("ui_viz_active_badge_pca"),
                downloadButton("download_pca_r_script_header", "Download R Script", class = "btn btn-outline-primary btn-sm", icon = icon("file-code")),
                span(class = "badge bg-light text-secondary border", "Metric Space")
              )
            )
          ),
          card_body(
            fluidRow(
              column(3, uiOutput("ui_pca_dataset_selector")),
              column(2, uiOutput("ui_pca_modality_selector")),
              column(3,
                     div(style = "margin-bottom: 8px;",
                       actionButton("btn_plot_pca", "Confirm & Plot PCA", class = "btn btn-primary btn-sm w-100", icon = icon("compass"), style = "font-weight: 600;")
                     ),
                     selectInput(
                       "sel_pca_cat", "Category Filter:",
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
              column(2, selectInput("sel_pca_panel", "PCA Panel:", choices = c("Both (Simulators & Metric Loadings)" = "both", "Simulator Ordination Only" = "methods", "Metric Loadings Only" = "loadings"), selected = "both")),
              column(2, numericInput("num_pca_top_metrics", "Top Metrics to Label:", value = 8, min = 3, max = 25, step = 1))
            ),
            fluidRow(
              column(12, style = "margin-top: 10px;",
                     downloadButton("download_pca_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2 mb-1"),
                     downloadButton("download_pca_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary me-2 mb-1"),
                     downloadButton("download_pca_r_script", "Download R Script", class = "btn btn-sm btn-outline-primary mb-1", icon = icon("file-code"))
              )
            ),
            fluidRow(
              column(6,
                     sliderInput("sld_pca_height", "PCA Height (px):", min = 500, max = 1600, value = 1100, step = 50)
              ),
              column(6,
                     sliderInput("sld_pca_width", "PCA Width (px):", min = 650, max = 1400, value = 950, step = 25)
              )
            ),
            div(
              style = "margin-bottom: 8px; font-size: 0.80rem; color: #64748B; background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 6px 12px;",
              icon("circle-info", class = "text-info me-1"),
              tags$b("Tip: "),
              "Arrows show the metrics that most strongly separate high-performing simulators from lower-performing ones."
            ),
            hr(),
            div(
              class = "bubble-scroll-container",
              style = "text-align: center; overflow-x: auto; padding: 10px;",
              uiOutput("ui_plot_metric_pca")
            )
          )
        )
      ),
      
      # Sub-panel 7: MDS Metric Space
      nav_panel(
        "7. MDS Metric Space",
        card(
          fill = FALSE,
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              span(icon("circle-nodes", class = "me-2 text-primary"), tags$b("MDS Distance Plot")),
              div(
                class = "d-flex align-items-center gap-2",
                uiOutput("ui_viz_active_badge_mds"),
                downloadButton("download_mds_r_script_header", "Download R Script", class = "btn btn-outline-primary btn-sm", icon = icon("file-code")),
                span(class = "badge bg-light text-secondary border", "Simulator Similarity")
              )
            )
          ),
          card_body(
            fluidRow(
              column(3, uiOutput("ui_mds_dataset_selector")),
              column(2, uiOutput("ui_mds_modality_selector")),
              column(3,
                     div(style = "margin-bottom: 8px;",
                       actionButton("btn_plot_mds", "Confirm & Plot MDS", class = "btn btn-primary btn-sm w-100", icon = icon("circle-nodes"), style = "font-weight: 600;")
                     ),
                     selectInput(
                       "sel_mds_cat", "Category Filter:",
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
              column(2, selectInput("sel_mds_by", "Compare By:", choices = c("By Simulators" = "simulators", "By Metric Summaries" = "summaries"), selected = "simulators")),
              column(2, sliderInput("sld_mds_height", "MDS Height (px):", min = 400, max = 1300, value = 620, step = 20))
            ),
            fluidRow(
              column(12, style = "margin-top: 10px;",
                     downloadButton("download_mds_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2 mb-1"),
                     downloadButton("download_mds_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary me-2 mb-1"),
                     downloadButton("download_mds_r_script", "Download R Script", class = "btn btn-sm btn-outline-primary mb-1", icon = icon("file-code"))
              )
            ),
            hr(),
            div(
              class = "bubble-scroll-container",
              style = "text-align: center; overflow-x: auto; padding: 10px;",
              uiOutput("ui_plot_metric_mds")
            )
          )
        )
      ),
      
      # Sub-panel 8: Cell Embeddings (t-SNE & UMAP)
      nav_panel(
        "8. Cell Embeddings (t-SNE & UMAP)",
        card(
          fill = FALSE,
          card_header(
            div(
              class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
              span(icon("project-diagram", class = "me-2 text-primary"), tags$b("Cell Projections (UMAP, t-SNE, PCA)")),
              div(
                class = "d-flex align-items-center gap-2",
                uiOutput("ui_viz_active_badge_emb"),
                downloadButton("download_emb_r_script_header", "Download R Script", class = "btn btn-outline-primary btn-sm", icon = icon("file-code")),
                span(class = "badge bg-light text-secondary border", "Cell Clusters")
              )
            )
          ),
          card_body(
            fluidRow(
              column(3,
                     uiOutput("ui_emb_dataset_selector"),
                     uiOutput("ui_emb_modality_selector"),
                     radioButtons(
                       "sel_emb_reduction", "1. Method:",
                       choices = c("UMAP" = "umap", "t-SNE" = "tsne", "PCA" = "pca"),
                       selected = "umap", inline = TRUE
                     ),
                     selectInput(
                       "sel_emb_layout", "2. Layout:",
                       choices = c(
                         "All Simulators Grid" = "facet",
                         "Side-by-Side (Reference vs Single Simulator)" = "side_by_side"
                       ),
                       selected = "facet"
                     ),
                     conditionalPanel(
                       condition = "input.sel_emb_layout == 'side_by_side'",
                       uiOutput("ui_emb_single_sim")
                     )
              ),
              column(3,
                     uiOutput("ui_emb_sim_picker"),
                     selectInput(
                       "sel_emb_color", "Color Cells By:",
                       choices = c(
                         "Cell Type / Group" = "cell_type",
                         "Cluster" = "cluster",
                         "Library Size (Depth)" = "library_size",
                         "Batch" = "batch",
                         "Dataset Source" = "dataset"
                       ),
                       selected = "cell_type"
                     )
              ),
              column(3,
                     sliderInput("sld_emb_pcs", "Number of PCs:", min = 5, max = 50, value = 20, step = 5),
                     conditionalPanel(
                       condition = "input.sel_emb_reduction == 'tsne'",
                       sliderInput("sld_emb_perp", "t-SNE Perplexity:", min = 5, max = 50, value = 15, step = 5)
                     ),
                     conditionalPanel(
                       condition = "input.sel_emb_reduction == 'umap'",
                       sliderInput("sld_emb_neighbors", "UMAP Neighbors:", min = 5, max = 50, value = 15, step = 5)
                     )
              ),
              column(3,
                     div(class = "mb-2",
                         actionButton("btn_plot_emb", "Compute & Render Embeddings", class = "btn btn-primary btn-sm w-100", icon = icon("project-diagram"), style = "font-weight: 600;")
                     ),
                     fluidRow(
                       column(6, sliderInput("sld_emb_pt_size", "Point Size:", min = 0.2, max = 3.0, value = 1.0, step = 0.1)),
                       column(6, sliderInput("sld_emb_alpha", "Alpha:", min = 0.2, max = 1.0, value = 0.8, step = 0.05))
                     ),
                     div(class = "mt-2",
                         downloadButton("download_emb_jpeg", "Download JPEG (600 DPI)", class = "btn btn-sm btn-primary me-2 mb-1"),
                         downloadButton("download_emb_pdf", "Download PDF", class = "btn btn-sm btn-outline-secondary me-2 mb-1"),
                         downloadButton("download_emb_r_script", "Download R Script", class = "btn btn-sm btn-outline-primary mb-1", icon = icon("file-code"))
                     )
              )
            ),
            hr(),
            div(
              style = "text-align: center; overflow-x: auto; padding: 10px;",
              uiOutput("ui_plot_cell_embeddings")
            ),
            hr(style = "margin: 20px 0; border-color: #CBD5E1;"),
            
            # Quantitative Metrics Box (Shown after user confirms and renders embeddings)
            conditionalPanel(
              condition = "input.btn_plot_emb > 0 || input.btn_plot_emb_inline > 0",
              div(
                class = "p-3", style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 8px;",
                div(
                  class = "d-flex justify-content-between align-items-center mb-2",
                  h6(icon("chart-line", class = "text-primary me-2"), tags$b("Cluster Separation and Quality Metrics"), style = "color: #0F172A; margin: 0; font-size: 0.92rem;"),
                  span(class = "badge bg-light text-secondary border", "Silhouette • ARI • Library Depth")
                ),
                p("Shows how well cell clusters and library sizes are preserved between real and simulated data:", style = "font-size: 0.82rem; color: #64748B; margin-bottom: 10px;"),
                DTOutput("table_emb_quality_metrics")
              )
            )
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 5: BENCHMARK RESULTS & EXPORT
  # ============================================================================
  nav_panel(
    "Benchmark Results & Export",
    # Benchmark Results Session Export (Exclusively .rds for reproducible analysis)
    div(
      class = "p-3 mb-3",
      style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
      div(
        class = "d-flex flex-wrap justify-content-between align-items-center gap-3",
        div(
          div(
            class = "d-flex align-items-center gap-2",
            icon("database", class = "text-primary fs-5"),
            tags$h5("Benchmark Results Session Export", style = "margin: 0; font-weight: 700; color: #0F172A; font-size: 1.05rem;"),
            span(class = "badge bg-primary-subtle text-primary border border-primary-subtle", "Complete R Data Object (.rds)")
          ),
          p("Download complete benchmarking session object (.rds) to preserve and reproduce all 62 evaluation metrics, rankings, and figures in another R session:",
            style = "margin: 2px 0 0 0; font-size: 0.82rem; color: #64748B;")
        ),
        div(
          downloadButton(
            "download_benchmark_data", "Download Benchmark Session (.rds)",
            class = "btn btn-primary", icon = icon("file-code"),
            style = "font-weight: 600; white-space: nowrap;"
          )
        )
      )
    ),
    
    # Master Results Table
    card(
      card_header(
        div(
          class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
          div(
            div(
              class = "d-flex align-items-center gap-2",
              icon("table", class = "text-primary"),
              tags$b("Master Benchmark Scores Table", style = "color: #1E293B;"),
              span(class = "badge bg-primary-subtle text-primary border border-primary-subtle", "Raw Metric Values & Normalized Scores")
            ),
            div(
              class = "text-muted small mt-1", style = "font-size: 0.79rem;",
              tags$span(tags$b("Score_Raw: "), "Raw statistical/computational values (e.g. KS distance, Wasserstein distance, runtime in seconds, peak RAM in MiB, Chi-Square values, correlation). "),
              tags$span(tags$b("Score: "), "Direction-aware normalized score [0, 1] where 1.0 represents best performance.")
            )
          ),
          div(
            class = "d-flex align-items-center gap-2",
            downloadButton("download_master_table_excel", "Download Table (.xlsx)", class = "btn btn-sm btn-outline-success", icon = icon("file-excel")),
            downloadButton("download_master_table_csv", "Download Table (.csv)", class = "btn btn-sm btn-outline-primary", icon = icon("file-csv")),
            span(class = "badge bg-light text-muted border", "Dropdown Filters Active")
          )
        )
      ),
      card_body(
        style = "padding: 12px;",
        DTOutput("table_master_export")
      )
    )
  ),
  
  # ============================================================================
  # TAB 7: HELP & GETTING STARTED (COMPREHENSIVELY UPDATED USER GUIDE)
  # ============================================================================
  nav_panel(
    "Help & Getting Started",
    fluidRow(
      column(
        12,
        card(
          card_body(
            # Introduction Hero Callout
            div(
              style = "background: #F8FAFC; border-left: 4px solid #1B4F72; border-radius: 6px; padding: 18px 22px; margin-bottom: 24px;",
              h4("Simulation Evaluation for Single-Cell and Multiomics Datasets", style = "font-weight: 800; color: #1B4F72; margin-top: 0; font-size: 1.18rem;"),
              p(
                strong("scSimEval"), " is a comprehensive, publication-grade benchmarking suite for single-cell RNA-seq (scRNA-seq), single-cell ATAC-seq (scATAC-seq), and single-cell multiomics (scRNA-seq + scATAC-seq) simulations. ",
                "It systematically evaluates how faithfully simulated datasets reproduce real experimental biology across ", strong("62 quantitative evaluation metrics organized into 8 canonical categories"), "."
              ),
              p(
                "This updated web application features multi-dataset benchmarking, interactive simulator leaderboards, customizable weight modeling, high-density bubble matrices, 8 diagnostic visualization suites, and reproducible R script generators for all plots.",
                style = "margin-bottom: 0; color: #475569; font-size: 0.88rem;"
              ),
              div(
                style = "margin-top: 14px; padding-top: 12px; border-top: 1px solid #CBD5E1; display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 10px;",
                span(
                  tags$i(class = "fa fa-book", style = "margin-right: 6px; color: #1B4F72;"),
                  strong("Online Documentation & Vignettes: "),
                  tags$a(
                    href = "https://kabilanbio.github.io/scSimEval",
                    target = "_blank",
                    rel = "noopener noreferrer",
                    "https://kabilanbio.github.io/scSimEval",
                    style = "color: #0284C7; text-decoration: underline; font-weight: 700; font-size: 0.90rem; margin-left: 4px;"
                  )
                ),
                tags$a(
                  href = "https://kabilanbio.github.io/scSimEval",
                  target = "_blank",
                  rel = "noopener noreferrer",
                  class = "btn btn-sm btn-primary",
                  style = "font-weight: 600; border-radius: 6px;",
                  tags$i(class = "fa fa-external-link-alt", style = "margin-right: 5px;"),
                  "Open Online Documentation"
                )
              )
            ),
            
            # ------------------------------------------------------------------
            # Section 1: Workflow Architecture
            # ------------------------------------------------------------------
            h4("1. Workflow Architecture & Pipeline", style = "font-weight: 700; color: #1B4F72; margin-top: 10px; font-size: 1.08rem;"),
            p("The benchmarking pipeline standardizes simulation evaluation across single and multi-dataset experiments:", style = "color: #475569; font-size: 0.88rem;"),
            
            div(
              style = "text-align: center; margin: 20px auto 16px auto; max-width: 980px;",
              tags$img(
                src = "scfigures/workflow_diagram.png",
                alt = "scSimEval Benchmarking Workflow Architecture",
                style = "width: 100%; max-width: 960px; height: auto; display: block; margin: 0 auto; border-radius: 8px; border: 1px solid #E2E8F0;"
              ),
              p(
                tags$b("Figure 1 | The scSimEval Benchmarking Workflow. "),
                "Step 1: Input real reference and simulated count matrices with metadata. Step 2: Calculate 62 evaluation metrics across 8 categories. Step 3: Standardize and consolidate all scores into a unified multi-dataset benchmark table. Step 4: Explore simulator rankings, bubble matrices, and high-resolution figures.",
                style = "font-size: 0.84rem; color: #64748B; margin-top: 10px; max-width: 960px; margin-left: auto; margin-right: auto;"
              )
            ),
            
            tags$div(
              class = "row g-3 my-2",
              tags$div(
                class = "col-md-6",
                tags$div(
                  style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("Step 1: Input Experimental Reference & Simulated Data"), style = "color: #1B4F72; font-size: 0.90rem;"),
                  p("Provide real reference count matrices and simulated datasets (genes x cells for scRNA-seq, peaks x cells for scATAC-seq), along with optional cell type annotations, batch labels, and computational resource records (runtime in seconds and peak RAM in MiB).", style = "font-size: 0.84rem; color: #475569; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-6",
                tags$div(
                  style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("Step 2: Calculate 62 Quantitative Metrics"), style = "color: #1B4F72; font-size: 0.90rem;"),
                  p("scSimEval computes 62 quantitative measures across 8 biological and computational categories, assessing statistical properties, zero-inflation, clustering concordance, biological fidelity, batch mixing, and scalability.", style = "font-size: 0.84rem; color: #475569; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-6",
                tags$div(
                  style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("Step 3: Direction-Aware Normalization & Consolidation"), style = "color: #1B4F72; font-size: 0.90rem;"),
                  p("The benchmark engine standardizes raw metric values into direction-aware normalized scores [0, 1] (where 1.0 is best) while preserving unnormalized raw values (Score_Raw) for full transparency.", style = "font-size: 0.84rem; color: #475569; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-6",
                tags$div(
                  style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("Step 4: Interactive Leaderboards, Visualizations & Standalone Scripts"), style = "color: #1B4F72; font-size: 0.90rem;"),
                  p("Explore cross-dataset bubble matrices, simulator ranking podiums, radar charts, distribution QC curves, PCA biplots, and MDS spaces. Download standalone R scripts for 100% reproducible figures outside Shiny.", style = "font-size: 0.84rem; color: #475569; margin-bottom: 0;")
                )
              )
            ),
            hr(style = "margin: 28px 0;"),
            
            # ------------------------------------------------------------------
            # Section 2: The Eight Evaluation Categories
            # ------------------------------------------------------------------
            h4("2. The Eight Evaluation Categories", style = "font-weight: 700; color: #1B4F72; font-size: 1.08rem;"),
            p("The 62 evaluation criteria in scSimEval are structured across eight foundational categories covering statistical, cellular, molecular, and computational dimensions:", style = "color: #475569; font-size: 0.88rem;"),
            
            tags$div(
              class = "table-responsive",
              tags$table(
                class = "table table-hover table-bordered table-sm",
                style = "font-size: 0.84rem;",
                tags$thead(
                  class = "table-light",
                  tags$tr(
                    tags$th(style = "width: 5%; text-align: center;", "#"),
                    tags$th(style = "width: 25%;", "Category"),
                    tags$th(style = "width: 15%; text-align: center;", "Focus"),
                    tags$th(style = "width: 55%;", "Representative Metrics & Scope")
                  )
                ),
                tags$tbody(
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "I"),
                    tags$td(tags$b("Distributional Properties")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-primary border", "Statistical")),
                    tags$td("Mean expression, variance, dispersion (CV, Fano factor), library size, detected feature rate, and distribution distances (Wasserstein, Kolmogorov-Smirnov statistic, Bhattacharyya distance, ECDF difference).")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "II"),
                    tags$td(tags$b("Correlations & Zero-Inflation")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-info border", "Co-expression")),
                    tags$td("Gene-gene Spearman/Pearson correlation networks, dropout rate curves, zero proportion concordance, and mean-variance trend fidelity.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "III"),
                    tags$td(tags$b("Cellular Structure & Concordance")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-success border", "Clustering")),
                    tags$td("Adjusted Rand Index (ARI), Normalized Mutual Information (NMI), Adjusted Mutual Information (AMI), Silhouette width, cluster purity, and cell-cell distance matrix concordance.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "IV"),
                    tags$td(tags$b("Batch Effects & Confounder Mixing")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-warning border", "Technical")),
                    tags$td("kBET (batch test), LISI (Local Inverse Simpson Index), Average Silhouette Width for Batch (ASW-batch), and principal component batch correlation.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "V"),
                    tags$td(tags$b("Biological Signal & Downstream Fidelity")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-danger border", "Phenotype")),
                    tags$td("Differentially Expressed Gene (DEG) concordance (Jaccard index, sensitivity, specificity, logFC correlation), marker gene ranking, and pathway enrichment overlap.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "VI"),
                    tags$td(tags$b("Trajectory & Lineage Dynamics")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-secondary border", "Pseudotime")),
                    tags$td("Pseudotime correlation with reference, trajectory topology preservation (tree height RMSE, graph alignment), and gene expression patterns along lineages.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "VII"),
                    tags$td(tags$b("Cross-Modal Coupling & Modularity")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-indigo border", "Multiomics")),
                    tags$td("Peak-to-gene linkage concordance, cross-modal cell pairing accuracy (FOSCTTM, Match@1), TF regulatory consistency, and multi-layer cluster concordance.")
                  ),
                  tags$tr(
                    tags$td(style = "text-align: center; font-weight: bold;", "VIII"),
                    tags$td(tags$b("Computational Scalability")),
                    tags$td(style = "text-align: center;", span(class = "badge bg-light text-dark border", "Efficiency")),
                    tags$td("Wall-clock runtime (seconds), peak RAM memory usage (MiB RAM), cell-scaling rate, and feature-scaling efficiency.")
                  )
                )
              )
            ),
            hr(style = "margin: 28px 0;"),
            
            # ------------------------------------------------------------------
            # Section 3: How Scores are Calculated and Normalized
            # ------------------------------------------------------------------
            h4("3. How Scores are Calculated and Normalized", style = "font-weight: 700; color: #1B4F72; font-size: 1.08rem;"),
            p("To fairly compare diverse metrics with differing physical units (e.g. KS distance, seconds, MiB, correlation), scSimEval standardizes all scores using a two-step normalization pipeline:", style = "color: #475569; font-size: 0.88rem;"),
            
            div(
              style = "text-align: center; margin: 20px auto 16px auto; max-width: 950px;",
              tags$img(
                src = "scfigures/score_normalization_workflow.png",
                alt = "Two-Step Score Normalization & Visual Mapping Pipeline",
                style = "width: 100%; max-width: 920px; height: auto; display: block; margin: 0 auto; border-radius: 8px; border: 1px solid #E2E8F0;"
              ),
              p(
                tags$b("Figure 2 | Score Normalization and Visual Mapping Workflow. "),
                "Metrics where lower values indicate better results (such as error, runtime, and memory) are inverted so that higher is always better. All values are scaled between 0.00 and 1.00 so simulators can be ranked and compared fairly. In the bubble matrix, bubble size reflects normalized fidelity, and top performers for each metric are highlighted with bold circular outlines.",
                style = "font-size: 0.84rem; color: #64748B; margin-top: 10px; max-width: 920px; margin-left: auto; margin-right: auto;"
              )
            ),
            
            tags$div(
              class = "row g-3 my-2",
              tags$div(
                class = "col-md-3",
                tags$div(
                  style = "background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("1. Direction Inversion"), style = "color: #1B4F72; font-size: 0.88rem;"),
                  p("Metrics where lower values are better (such as distances, errors, runtime in seconds, or peak RAM in MiB) are inverted:", style = "font-size: 0.80rem; color: #475569;"),
                  tags$p(tags$code("Inverted = 1 - (Value - Min) / (Max - Min)"), style = "text-align: center; font-weight: bold; font-size: 0.80rem;"),
                  p("This guarantees that higher numbers always signify higher simulation fidelity.", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-3",
                tags$div(
                  style = "background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("2. Scale Between 0.0 and 1.0"), style = "color: #1B4F72; font-size: 0.88rem;"),
                  p("All metrics are rescaled across simulation methods onto a strictly bounded [0.0, 1.0] scale:", style = "font-size: 0.80rem; color: #475569;"),
                  tags$p(tags$code("Score = (Value - Min) / (Max - Min)"), style = "text-align: center; font-weight: bold; font-size: 0.80rem;"),
                  p("0.0 represents the lowest performing simulator, while 1.0 represents top fidelity.", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-3",
                tags$div(
                  style = "background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("3. Composite Overall Score"), style = "color: #1B4F72; font-size: 0.88rem;"),
                  p("Overall simulator fidelity is calculated by taking the average (or user-weighted sum) of normalized scores:", style = "font-size: 0.80rem; color: #475569;"),
                  tags$p(tags$code("Composite = sum(w_i * Score_i)"), style = "text-align: center; font-weight: bold; font-size: 0.80rem;"),
                  p("Displayed as a score from 0.0 to 1.0 and as a percentage on the Leaderboard.", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 0;")
                )
              ),
              tags$div(
                class = "col-md-3",
                tags$div(
                  style = "background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 8px; padding: 14px 16px; height: 100%;",
                  h6(tags$b("4. Dual Value Reporting"), style = "color: #1B4F72; font-size: 0.88rem;"),
                  p("The Master Benchmark Scores Table provides both representation forms:", style = "font-size: 0.80rem; color: #475569;"),
                  tags$ul(
                    style = "font-size: 0.78rem; color: #475569; padding-left: 14px; margin-bottom: 0;",
                    tags$li(tags$b("Score_Raw: "), "Raw statistical metric (seconds, MiB, KS statistic)."),
                    tags$li(tags$b("Score: "), "Normalized score [0, 1] used for rankings & bubble sizes.")
                  )
                )
              )
            ),
            hr(style = "margin: 28px 0;"),
            
            # ------------------------------------------------------------------
            # Section 4: Comprehensive Guide to the 8 App Tabs
            # ------------------------------------------------------------------
            h4("4. Comprehensive Guide to the 8 App Tabs", style = "font-weight: 700; color: #1B4F72; font-size: 1.08rem;"),
            p("Navigate through the application tabs to perform simulation evaluations and explore benchmarks:", style = "color: #475569; font-size: 0.88rem;"),
            
            tags$div(
              class = "accordion my-3", id = "accordionGuide",
              
              # Tab 1: Home
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingOne",
                  tags$button(
                    class = "accordion-button", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseOne", "aria-expanded" = "true", "aria-controls" = "collapseOne",
                    tags$b("Tab 1: Home — Framework Overview & Operational Modes")
                  )
                ),
                tags$div(
                  id = "collapseOne", class = "accordion-collapse collapse show", "aria-labelledby" = "headingOne",
                  tags$div(
                    class = "accordion-body",
                    p("Presents an executive overview of scSimEval, metric counts across all 8 biological and computational categories, and an interactive operational mode tour showcasing Pre-computed Benchmarks, Single-Cell RNA-seq Evaluation, Multi-Omics Evaluation, and Pre-computed Table Upload.", style = "font-size: 0.84rem; color: #475569; margin: 0;")
                  )
                )
              ),
              
              # Tab 2: Data Hub
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingTwo",
                  tags$button(
                    class = "accordion-button collapsed", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseTwo", "aria-expanded" = "false", "aria-controls" = "collapseTwo",
                    tags$b("Tab 2: Data Hub — Multi-Dataset Ingestion & Batch HPC Studio")
                  )
                ),
                tags$div(
                  id = "collapseTwo", class = "accordion-collapse collapse", "aria-labelledby" = "headingTwo",
                  tags$div(
                    class = "accordion-body",
                    tags$ul(
                      style = "font-size: 0.84rem; color: #475569; padding-left: 16px; margin: 0;",
                      tags$li(tags$b("Pre-Computed Benchmark (Demo Data): "), "Load a precomputed 6-simulator multi-dataset benchmark across Dataset 1 and Dataset 2 for instant exploration."),
                      tags$li(tags$b("Single-Cell RNA-seq Evaluation: "), "Upload real reference matrices and simulated datasets (.rds, .csv, .tsv, .txt) with optional cell type annotations and batch labels to calculate all 62 evaluation metrics."),
                      tags$li(tags$b("Multi-Omics (RNA + ATAC) Evaluation: "), "Benchmark paired single-cell RNA-seq and ATAC-seq co-assays evaluating cross-modal coupling, peak-to-gene linkages, and feature alignment."),
                      tags$li(tags$b("Upload Pre-Computed Benchmark Table: "), "Directly ingest previously calculated benchmark CSV/RDS tables for instant visualization and leaderboard generation."),
                      tags$li(tags$b("Multi-Dataset Benchmarking: "), "Upload or evaluate across multiple experimental conditions (Dataset 1, Dataset 2, ..., Dataset N) with dynamic simulator renaming."),
                      tags$li(tags$b("Dataset Summary Table: "), "Inspect comprehensive QC statistics across all datasets (cell count, gene count, sparsity %, sequencing depth, dropout rates)."),
                      tags$li(tags$b("HPC Cluster Batch Script Generator: "), "Export automated SLURM, PBS, or Bash cluster execution scripts to run computationally intensive evaluations offline.")
                    )
                  )
                )
              ),
              
              # Tab 3: Leaderboard
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingThree",
                  tags$button(
                    class = "accordion-button collapsed", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseThree", "aria-expanded" = "false", "aria-controls" = "collapseThree",
                    tags$b("Tab 3: Leaderboard & Simulator Rankings — Composite Podiums & Weight Customization")
                  )
                ),
                tags$div(
                  id = "collapseThree", class = "accordion-collapse collapse", "aria-labelledby" = "headingThree",
                  tags$div(
                    class = "accordion-body",
                    tags$ul(
                      style = "font-size: 0.84rem; color: #475569; padding-left: 16px; margin: 0;",
                      tags$li(tags$b("Simulator Ranking Cards & Podiums: "), "Visual podium ranking simulators by overall composite fidelity across all 8 evaluation categories."),
                      tags$li(tags$b("Category Score Breakdown: "), "Detailed breakdown of performance per simulator across each individual category (Distributional Properties, Clustering, Biological Signal, etc.)."),
                      tags$li(tags$b("Radar / Star Composite Plots: "), "Multi-axis radial plots highlighting simulator trade-offs, strengths, and weaknesses."),
                      tags$li(tags$b("Custom Category Weight Sliders: "), "Interactively adjust the relative importance of each category to recalculate composite rankings according to specific biological priorities.")
                    )
                  )
                )
              ),
              
              # Tab 4: Bubble Matrix
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingFour",
                  tags$button(
                    class = "accordion-button collapsed", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseFour", "aria-expanded" = "false", "aria-controls" = "collapseFour",
                    tags$b("Tab 4: Comparative Bubble Matrix — High-Density Cross-Dataset Benchmark Matrix")
                  )
                ),
                tags$div(
                  id = "collapseFour", class = "accordion-collapse collapse", "aria-labelledby" = "headingFour",
                  tags$div(
                    class = "accordion-body",
                    tags$ul(
                      style = "font-size: 0.84rem; color: #475569; padding-left: 16px; margin: 0;",
                      tags$li(tags$b("Multi-Dataset Selector Bar: "), "Switch seamlessly between All Datasets (Cross-Dataset Average) or inspect individual datasets (Dataset 1, Dataset 2, ..., Dataset N) separately."),
                      tags$li(tags$b("Minimalist vs Publication Matrix: "), "Toggle between a sleek, dot-matrix view and a fully annotated publication matrix."),
                      tags$li(tags$b("Category & Modality Filtering: "), "Filter metrics by modality (All Categories, RNA-Specific, ATAC-Specific, or Cross-Modal Coupling)."),
                      tags$li(tags$b("High-Resolution Exports: "), "Export publication figures in 600 DPI JPEG or vectorized PDF."),
                      tags$li(tags$b("Download Standalone R Script: "), "Generate a 100% self-contained, end-to-end analytical R script without readymade dput dumps, featuring a dedicated FILE PATH CONFIGURATION & HPC CLUSTER RUN GUIDE to run on local workstations, cloud instances, or SLURM/PBS clusters.")
                    )
                  )
                )
              ),
              
              # Tab 5: Visualizations
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingFive",
                  tags$button(
                    class = "accordion-button collapsed", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseFive", "aria-expanded" = "false", "aria-controls" = "collapseFive",
                    tags$b("Tab 5: Visualizations — 8 Dedicated Diagnostic Panels & Reproducible R Scripts")
                  )
                ),
                tags$div(
                  id = "collapseFive", class = "accordion-collapse collapse", "aria-labelledby" = "headingFive",
                  tags$div(
                    class = "accordion-body",
                    tags$p("Explore eight dedicated diagnostic subtabs equipped with global multi-dataset pill selectors, modality filtering (scRNA-seq, scATAC-seq, Cross-Modal Coupling), and standalone reproducible R script generators:", style = "font-size: 0.84rem; color: #475569; margin-bottom: 6px;"),
                    tags$ol(
                      style = "font-size: 0.84rem; color: #475569; padding-left: 18px; margin: 0;",
                      tags$li(tags$b("Evaluation Summary: "), "Horizontal bar charts ranking simulators across categories and composite fidelity with dynamic raw vs. normalized score scaling."),
                      tags$li(tags$b("Distribution QC: "), "Empirical cumulative distribution curves (ECDF), KS distance curves, and gene/cell mean-variance trends comparing reference vs simulations with modality isolation."),
                      tags$li(tags$b("Scalability Benchmark: "), "Runtime in seconds vs Peak RAM in MiB, throughput (cells/sec), and resource cost index."),
                      tags$li(tags$b("Metric Boxplots & Barplots: "), "Score distribution across categories with Raw vs Normalized view toggles."),
                      tags$li(tags$b("Metric Heatmap: "), "Unsupervised hierarchical clustering of simulators and metrics with dynamic score scaling and custom dimensions."),
                      tags$li(tags$b("PCA Ordination: "), "Principal component biplots projecting simulators into reduced metric space."),
                      tags$li(tags$b("MDS Metric Space: "), "Multidimensional scaling mapping simulator algorithmic proximity."),
                      tags$li(tags$b("Cell Embeddings: "), "Dimension reduction (UMAP, t-SNE, PCA) comparing reference cells directly to simulated cells with modality isolation and quantitative fidelity metrics.")
                    )
                  )
                )
              ),
              
              # Tab 6: Benchmark Results & Export
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingSix",
                  tags$button(
                    class = "accordion-button collapsed", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseSix", "aria-expanded" = "false", "aria-controls" = "collapseSix",
                    tags$b("Tab 6: Benchmark Results & Export — Master Scores & Full Data Archives")
                  )
                ),
                tags$div(
                  id = "collapseSix", class = "accordion-collapse collapse", "aria-labelledby" = "headingSix",
                  tags$div(
                    class = "accordion-body",
                    tags$ul(
                      style = "font-size: 0.84rem; color: #475569; padding-left: 16px; margin: 0;",
                      tags$li(tags$b("Master Benchmark Scores Table: "), "Interactive, searchable, and sortable table with instant column-level dropdown filtering for Dataset, Method, and Category."),
                      tags$li(tags$b("Raw & Normalized Values: "), "Explicitly reports Score_Raw (raw physical/statistical metric measurements, runtime in seconds, peak RAM in MiB) alongside Score (normalized fidelity score [0, 1])."),
                      tags$li(tags$b("Multi-Format Exports: "), "Export complete benchmark results as an Excel workbook (.xlsx with individual sheets), CSV table (.csv), Tab-delimited file (.tsv), R object (.rds), or complete ZIP archive (.zip with README).")
                    )
                  )
                )
              ),
              
              # Tab 7: Help & Getting Started
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingSeven",
                  tags$button(
                    class = "accordion-button collapsed", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseSeven", "aria-expanded" = "false", "aria-controls" = "collapseSeven",
                    tags$b("Tab 7: Help & Getting Started — Comprehensive User Guide & Formulas")
                  )
                ),
                tags$div(
                  id = "collapseSeven", class = "accordion-collapse collapse", "aria-labelledby" = "headingSeven",
                  tags$div(
                    class = "accordion-body",
                    p("Step-by-step documentation, metric glossary, mathematical normalization formulas, hardware guidelines, and input file templates.", style = "font-size: 0.84rem; color: #475569; margin: 0;")
                  )
                )
              ),
              
              # Tab 8: Contact
              tags$div(
                class = "accordion-item",
                tags$h2(
                  class = "accordion-header", id = "headingEight",
                  tags$button(
                    class = "accordion-button collapsed", type = "button", "data-bs-toggle" = "collapse", "data-bs-target" = "#collapseEight", "aria-expanded" = "false", "aria-controls" = "collapseEight",
                    tags$b("Tab 8: Contact — Laboratory & Citation Information")
                  )
                ),
                tags$div(
                  id = "collapseEight", class = "accordion-collapse collapse", "aria-labelledby" = "headingEight",
                  tags$div(
                    class = "accordion-body",
                    p("Information about the development team at Division of Agricultural Bioinformatics, ICAR-IASRI, maintainer email contacts, GitHub issues link, and recommended citations.", style = "font-size: 0.84rem; color: #475569; margin: 0;")
                  )
                )
              )
            ),
            
            # ------------------------------------------------------------------
            # Section 5: Input Data Formats & Metadata Templates
            # ------------------------------------------------------------------
            div(
              style = "margin-top: 24px; padding: 22px 24px; background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
              div(
                class = "d-flex justify-content-between align-items-center mb-2",
                h5(tags$b("Input Data Formats & Metadata Templates"), style = "color: #1E3A8A; margin: 0; font-size: 1.05rem;"),
                span(class = "badge bg-light text-primary border", "Required & Optional Inputs")
              ),
              p("scSimEval accepts standard single-cell matrix and annotation files in R object (.rds) or comma-/tab-delimited (.csv / .tsv / .txt) text formats:", style = "color: #475569; font-size: 0.85rem; margin-bottom: 16px;"),
              
              fluidRow(
                column(
                  4,
                  div(
                    style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 14px; height: 100%;",
                    h6(tags$b(icon("table-cells", class = "me-1 text-primary"), " 1. Count Matrices"), style = "color: #0F172A; margin-bottom: 6px; font-size: 0.88rem;"),
                    p(tags$b("Orientation: "), "Genes/Peaks in rows, Cells in columns.", style = "font-size: 0.79rem; color: #475569; margin-bottom: 6px;"),
                    p(tags$b("Formats: "), ".rds (dgCMatrix, matrix, SingleCellExperiment, Seurat) or .csv / .tsv with gene names in first column and cell barcodes in header.", style = "font-size: 0.77rem; color: #64748B; margin-bottom: 0;")
                  )
                ),
                column(
                  4,
                  div(
                    style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 14px; height: 100%; display: flex; flex-direction: column; justify-content: space-between;",
                    div(
                      h6(tags$b(icon("layer-group", class = "me-1 text-primary"), " 2. Cell Types (Optional)"), style = "color: #0F172A; margin-bottom: 6px; font-size: 0.88rem;"),
                      p("Used for clustering fidelity (ARI, NMI, Silhouette), marker gene DEGs, and embedding colors.", style = "font-size: 0.77rem; color: #64748B; margin-bottom: 8px;"),
                      tags$pre("Cell_ID,Cell_Type
Cell_01,TypeA
Cell_02,TypeA
Cell_03,TypeB
...", style = "font-size: 0.73rem; background: #FFFFFF; border: 1px solid #CBD5E1; padding: 6px 8px; border-radius: 4px; margin-bottom: 10px;")
                    ),
                    downloadButton("dl_help_sample_celltypes_csv", "Download Example Cell Types (.csv)", class = "btn btn-sm btn-outline-primary w-100", icon = icon("download"))
                  )
                ),
                column(
                  4,
                  div(
                    style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 14px; height: 100%; display: flex; flex-direction: column; justify-content: space-between;",
                    div(
                      h6(tags$b(icon("boxes-stacked", class = "me-1 text-teal"), " 3. Batch Labels (Optional)"), style = "color: #0F172A; margin-bottom: 6px; font-size: 0.88rem;"),
                      p("Used for Category IV batch mixing (kBET, CMS, LISI, variance decomposition).", style = "font-size: 0.77rem; color: #64748B; margin-bottom: 8px;"),
                      tags$pre("Cell_ID,Batch
Cell_01,Batch1
Cell_02,Batch2
Cell_03,Batch1
...", style = "font-size: 0.73rem; background: #FFFFFF; border: 1px solid #CBD5E1; padding: 6px 8px; border-radius: 4px; margin-bottom: 10px;")
                    ),
                    downloadButton("dl_help_sample_batch_csv", "Download Example Batch (.csv)", class = "btn btn-sm btn-outline-success w-100", icon = icon("download"))
                  )
                )
              )
            ),
            
            # Section 5: Online Documentation Box
            div(
              style = "margin-top: 24px; padding: 18px 20px; background: #F0F9FF; border: 1px solid #BAE6FD; border-radius: 8px;",
              div(
                class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
                div(
                  h6(tags$b("Online Documentation, Tutorials & Guides"), style = "color: #0369A1; margin: 0; font-size: 0.95rem;"),
                  p("Access comprehensive Vignettes, tutorials, and function reference at GitHub Pages:", style = "font-size: 0.82rem; color: #0284C7; margin: 2px 0 0 0;")
                ),
                tags$a(
                  href = "https://kabilanbio.github.io/scSimEval",
                  target = "_blank",
                  rel = "noopener noreferrer",
                  class = "btn btn-sm btn-outline-primary",
                  style = "font-weight: 600;",
                  tags$i(class = "fa fa-external-link-alt me-1"), "https://kabilanbio.github.io/scSimEval"
                )
              )
            ),
            
            # ------------------------------------------------------------------
            # Section 6: End-to-End Analytical R Scripts & HPC Cluster Execution Guide
            # ------------------------------------------------------------------
            div(
              style = "margin-top: 24px; padding: 22px 24px; background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
              div(
                class = "d-flex justify-content-between align-items-center mb-2 flex-wrap gap-2",
                h5(tags$b("End-to-End Analytical R Scripts & HPC Cluster Run Guide"), style = "color: #1E3A8A; margin: 0; font-size: 1.05rem;"),
                span(class = "badge bg-primary-subtle text-primary border border-primary-subtle", "100% Reproducible • Cloud & HPC Ready")
              ),
              p(
                "Every visualization panel in scSimEval includes a ", tags$b("Download Standalone R Script"), 
                " button. Unlike conventional web applications that dump hardcoded static numbers (e.g. ", tags$code("dput"), 
                "), scSimEval generates ", tags$b("complete, self-contained end-to-end analytical workflows"), 
                " that load raw input files from disk, execute full statistical calculations, and export publication-ready graphics.",
                style = "color: #475569; font-size: 0.85rem; margin-bottom: 16px;"
              ),
              
              fluidRow(
                column(
                  4,
                  div(
                    style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 14px; height: 100%;",
                    h6(tags$b(icon("code", class = "me-1 text-primary"), " 1. Procedural & dput-Free"), style = "color: #0F172A; margin-bottom: 6px; font-size: 0.88rem;"),
                    p("Scripts contain zero hardcoded data matrices. Clean, human-readable code that can be tracked in Git, peer-reviewed, and customized for manuscripts.", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 0;")
                  )
                ),
                column(
                  4,
                  div(
                    style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 14px; height: 100%;",
                    h6(tags$b(icon("sliders", class = "me-1 text-teal"), " 2. Dedicated Configuration Block"), style = "color: #0F172A; margin-bottom: 6px; font-size: 0.88rem;"),
                    p("Every script concludes with a dedicated FILE PATH CONFIGURATION block where non-coding users simply update input file paths before execution.", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 0;")
                  )
                ),
                column(
                  4,
                  div(
                    style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 14px; height: 100%;",
                    h6(tags$b(icon("server", class = "me-1 text-indigo"), " 3. HPC & Cloud Compatible"), style = "color: #0F172A; margin-bottom: 6px; font-size: 0.88rem;"),
                    p("Headless execution via Rscript. Run directly on SLURM, PBS, LSF, AWS EC2, GCP, or local machines without needing a web browser.", style = "font-size: 0.78rem; color: #64748B; margin-bottom: 0;")
                  )
                )
              ),
              
              div(
                style = "margin-top: 16px; background: #0F172A; border-radius: 6px; padding: 14px 18px; color: #E2E8F0;",
                div(
                  class = "d-flex justify-content-between align-items-center mb-2",
                  tags$span(tags$b(icon("terminal", class = "me-1 text-warning"), " Standardized Configuration Block Template (Bottom of Each Downloaded Script):"), style = "font-size: 0.80rem; font-weight: 600; color: #38BDF8;"),
                  tags$span(class = "badge bg-secondary", "Edit paths & run")
                ),
                tags$pre(
                  style = "font-size: 0.74rem; background: transparent; color: #F1F5F9; border: none; padding: 0; margin: 0; line-height: 1.5;",
                  "# ==============================================================================\n# FILE PATH CONFIGURATION & HPC CLUSTER RUN GUIDE\n# ==============================================================================\nBENCHMARK_RESULTS_FILE <- \"path/to/scSimEval_benchmark_results.csv\"  # Master results table\nREF_DATA_FILE          <- \"path/to/reference_counts.rds\"            # Empirical count matrix\nSIM_DATA_FILE          <- \"path/to/simulated_counts.rds\"            # Simulator count matrix\nCELL_TYPES_FILE        <- \"path/to/cell_types.csv\"                  # Optional cell annotations\nBATCH_FILE             <- \"path/to/batch_labels.csv\"                # Optional batch labels\nOUTPUT_DIR             <- \"./scSimEval_publication_figures\"         # Destination directory\n\n# Terminal execution on local workstation, cloud VM, or HPC node:\n#   $ Rscript scSimEval_workflow.R\n#\n# SLURM cluster submission:\n#   $ sbatch --mem=32G --cpus-per-task=4 --wrap=\"Rscript scSimEval_workflow.R\""
                )
              )
            ),
            
            # ------------------------------------------------------------------
            # Section 7: Hardware & Memory Guidelines
            # ------------------------------------------------------------------
            div(
              style = "margin-top: 22px; padding: 20px 22px; background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
              div(
                class = "d-flex justify-content-between align-items-center mb-2",
                h5(tags$b("Hardware & Memory Guidelines"), style = "color: #1E3A8A; margin: 0; font-size: 1.00rem;"),
                span(class = "badge bg-light text-primary border", "CPU & RAM Sizing")
              ),
              p(
                HTML("<b>Standard Analysis (Exploratory / Demo):</b> For datasets up to ~2,000 cells &times; 3,000 features, a standard workstation or cloud instance with <b>2–4 CPU cores and 8–16 GiB RAM</b> delivers rapid execution (<1 min) with minimal memory footprint."),
                style = "color: #334155; font-size: 0.85rem; margin-bottom: 6px; line-height: 1.5;"
              ),
              p(
                HTML("<b>Large-Scale & Multiomics Workflows:</b> For large single-cell datasets (>5,000 cells) or joint multiomics assays (scRNA-seq + scATAC-seq with multiple simulators), <b>4–8 CPU cores and 16–32 GiB RAM</b> are recommended to accommodate parallel metric evaluations and latent space projections."),
                style = "color: #334155; font-size: 0.85rem; margin-bottom: 6px; line-height: 1.5;"
              ),
              p(
                HTML("<b>Atlas-Scale Datasets (HPC):</b> For atlas-scale cohorts (>50,000 cells), allocating <b>8+ CPU cores and 32–64+ GiB RAM</b> on an HPC cluster node or dedicated workstation ensures optimal performance using scSimEval's sparse matrix operations and multi-core parallelization."),
                style = "color: #334155; font-size: 0.85rem; margin-bottom: 12px; line-height: 1.5;"
              ),
              tags$div(
                style = "margin-top: 10px; background: #F8FAFC; border-radius: 6px; padding: 10px 14px; border: 1px solid #E2E8F0;",
                span(tags$b("Launch locally in R: "), style = "font-size: 0.82rem; color: #334155;"),
                tags$code("remotes::install_github('kabilanbio/scSimEval'); scSimEval::launch_scSimEval_app()", style = "font-size: 0.82rem; background: #FFFFFF; padding: 2px 6px; border: 1px solid #CBD5E1; border-radius: 4px; margin-left: 6px;")
              ),
              div(
                style = "margin-top: 10px; padding: 10px 14px; background: #F8FAFC; border-left: 3px solid #6366F1; border-radius: 4px; font-size: 0.82rem; color: #475569;",
                tags$b("Memory Optimization Note: "),
                "scSimEval utilizes sparse matrix arithmetic (dgCMatrix) and fast low-rank linear algebra to minimize peak memory consumption across all 62 evaluation metrics."
              )
            )
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # TAB 8: LOGS (LIVE EXECUTION LOGS & SYSTEM DIAGNOSTICS)
  # ============================================================================
  nav_panel(
    "Logs",
    div(
      style = "max-width: 1400px; margin: 0 auto; padding: 10px 0 30px 0;",
      
      # 1. Hero Header Banner
      div(
        class = "hero-box",
        style = "background: linear-gradient(135deg, #0F172A 0%, #1E293B 100%); border-left: 5px solid #38BDF8; margin-bottom: 20px; padding: 22px 28px; border-radius: 10px;",
        div(
          class = "d-flex justify-content-between align-items-center flex-wrap gap-3",
          div(
            h3(
              tags$span(icon("terminal", class = "me-2 text-info")),
              "Live Execution Logs & System Diagnostics",
              style = "margin: 0; font-weight: 800; color: #FFFFFF; font-size: 1.45rem;"
            ),
            p(
              "Real-time diagnostic stream capturing data loading, metric computation across 62 evaluation criteria, visualization rendering, memory profiling, and live ERRORS, WARNINGS, and NOTES.",
              style = "margin: 6px 0 0 0; color: #94A3B8; font-size: 0.88rem;"
            )
          ),
          div(
            class = "d-flex align-items-center gap-2",
            span(
              class = "badge",
              style = "background: rgba(16, 185, 129, 0.2); color: #34D399; border: 1px solid rgba(16, 185, 129, 0.4); padding: 7px 14px; font-size: 0.82rem; font-weight: 700; border-radius: 20px;",
              icon("circle", class = "me-1 text-success", style = "font-size: 0.65rem; animation: pulseLive 1.5s infinite;"),
              "LIVE MONITOR ACTIVE"
            ),
            span(
              class = "badge",
              style = "background: rgba(56, 189, 248, 0.15); color: #7DD3FC; border: 1px solid rgba(56, 189, 248, 0.3); padding: 7px 12px; font-size: 0.82rem; font-weight: 600; border-radius: 20px;",
              "scSimEval v0.99.9"
            )
          )
        )
      ),
      
      # 2. KPI Summary Badges Row
      fluidRow(
        column(
          width = 3,
          div(
            class = "stat-card",
            style = "border-left-color: #38BDF8; background: #FFFFFF;",
            div(class = "d-flex justify-content-between align-items-center",
              div(
                div(class = "stat-number text-dark", uiOutput("log_stat_total", inline = TRUE)),
                div(class = "stat-label", "Total Log Events")
              ),
              icon("list-check", class = "text-info", style = "font-size: 2.2rem; opacity: 0.7;")
            )
          )
        ),
        column(
          width = 3,
          div(
            class = "stat-card",
            style = "border-left-color: #EF4444; background: #FFFFFF;",
            div(class = "d-flex justify-content-between align-items-center",
              div(
                div(class = "stat-number text-danger", uiOutput("log_stat_error", inline = TRUE)),
                div(class = "stat-label", "Critical Errors")
              ),
              icon("triangle-exclamation", class = "text-danger", style = "font-size: 2.2rem; opacity: 0.7;")
            )
          )
        ),
        column(
          width = 3,
          div(
            class = "stat-card",
            style = "border-left-color: #F59E0B; background: #FFFFFF;",
            div(class = "d-flex justify-content-between align-items-center",
              div(
                div(class = "stat-number text-warning", uiOutput("log_stat_warning", inline = TRUE)),
                div(class = "stat-label", "Warnings & Alerts")
              ),
              icon("circle-exclamation", class = "text-warning", style = "font-size: 2.2rem; opacity: 0.7;")
            )
          )
        ),
        column(
          width = 3,
          div(
            class = "stat-card",
            style = "border-left-color: #10B981; background: #FFFFFF;",
            div(class = "d-flex justify-content-between align-items-center",
              div(
                div(class = "stat-number text-success", uiOutput("log_stat_note_info", inline = TRUE)),
                div(class = "stat-label", "Notes, Info & Success")
              ),
              icon("circle-check", class = "text-success", style = "font-size: 2.2rem; opacity: 0.7;")
            )
          )
        )
      ),
      
      # 3. Interactive Toolbar Card
      card(
        class = "mb-3",
        style = "border: 1px solid #E2E8F0; border-radius: 8px;",
        card_body(
          style = "padding: 14px 18px;",
          div(
            class = "d-flex justify-content-between align-items-center flex-wrap gap-3",
            
            # Left: Level Filters
            div(
              class = "d-flex align-items-center flex-wrap gap-2",
              tags$span(tags$strong("Filter:"), style = "font-size: 0.85rem; color: #475569; margin-right: 4px;"),
              radioButtons(
                inputId = "log_filter_level",
                label = NULL,
                choices = c("All" = "ALL", "Errors" = "ERROR", "Warnings" = "WARNING", "Notes" = "NOTE", "Info" = "INFO", "Success" = "SUCCESS"),
                selected = "ALL",
                inline = TRUE
              )
            ),
            
            # Middle: Search Box
            div(
              style = "min-width: 260px; flex-grow: 1; max-width: 380px;",
              textInput(
                inputId = "log_search_query",
                label = NULL,
                placeholder = "Search log messages, components, errors...",
                width = "100%"
              )
            ),
            
            # Right: Controls & Export
            div(
              class = "d-flex align-items-center gap-2",
              checkboxInput("log_autoscroll", "Auto-scroll to bottom", value = TRUE, width = "auto"),
              actionButton(
                inputId = "btn_refresh_logs",
                label = "Refresh",
                icon = icon("arrows-rotate"),
                class = "btn btn-outline-secondary btn-sm",
                style = "font-weight: 600;"
              ),
              actionButton(
                inputId = "btn_clear_logs",
                label = "Clear Logs",
                icon = icon("trash-can"),
                class = "btn btn-outline-danger btn-sm",
                style = "font-weight: 600;"
              ),
              downloadButton(
                outputId = "btn_download_logs",
                label = "Download Logs (.txt)",
                class = "btn btn-primary btn-sm",
                icon = icon("download"),
                style = "font-weight: 600;"
              )
            )
          )
        )
      ),
      
      # 4. Live Terminal / Console Window
      div(
        class = "log-terminal-window",
        div(
          class = "log-terminal-header",
          div(
            class = "d-flex align-items-center gap-2",
            div(class = "terminal-dot terminal-dot-red"),
            div(class = "terminal-dot terminal-dot-yellow"),
            div(class = "terminal-dot terminal-dot-green"),
            span("scSimEval Diagnostics Console", style = "margin-left: 8px; font-weight: 600; letter-spacing: 0.3px; color: #E2E8F0; font-size: 0.82rem;"),
            span("● live output", style = "font-size: 0.72rem; color: #10B981; font-weight: 600; margin-left: 4px;")
          ),
          div(
            class = "d-flex align-items-center gap-3",
            span(uiOutput("log_console_status_line", inline = TRUE), style = "font-size: 0.78rem; color: #94A3B8;"),
            tags$button(
              type = "button",
              class = "btn btn-sm text-light p-0",
              style = "opacity: 0.75; font-size: 0.8rem; background: none; border: none; cursor: pointer;",
              onclick = "navigator.clipboard.writeText(document.getElementById('sc_log_body').innerText); alert('Logs copied to clipboard!');",
              title = "Copy raw console logs to clipboard",
              icon("copy", class = "me-1"), "Copy"
            )
          )
        ),
        div(
          id = "sc_log_body",
          class = "log-terminal-body",
          uiOutput("ui_log_console")
        )
      ),
      
      # 5. Explanatory Footer Callout
      div(
        style = "margin-top: 18px; padding: 12px 18px; background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 8px; font-size: 0.82rem; color: #64748B;",
        div(
          class = "d-flex align-items-center gap-3 flex-wrap",
          div(strong(icon("info-circle", class = "me-1 text-primary"), "Live Logging System Architecture: ")),
          div(span(class = "badge", style = "background: #7F1D1D; color: #FCA5A5;", "ERROR"), " Highlights computation failures, dimension mismatches, missing required packages, or execution interrupts."),
          div(span(class = "badge", style = "background: #78350F; color: #FDE68A;", "WARNING"), " Tracks metric boundary alerts, missing optional cell types/batches, or fallback approximations."),
          div(span(class = "badge", style = "background: #4C1D95; color: #DDD6FE;", "NOTE"), " Documents data preprocessing decisions, sparse conversions, and category metric count unlock events."),
          div(span(class = "badge", style = "background: #0C4A6E; color: #7DD3FC;", "INFO"), " Continuous timeline of user clicks, plot rendering triggers, dataset switches, and export generations.")
        )
      )
    )
  ),
    # ============================================================================
  # TAB 9: CONTACT (SINGLE ICAR LOGO, PRINCIPAL SCIENTIST, KABILAN SAKTHIVEL)
  # ============================================================================
  nav_panel(
    "Contact",
    fluidRow(
      column(
        12,
        card(
          card_body(
            # Institute Header Section (Single ICAR Logo)
            div(
              style = "margin-bottom: 22px; padding: 18px 24px; background: #FFFFFF; border-radius: 10px; border: 1px solid #E2E8F0; box-shadow: 0 1px 4px rgba(0,0,0,0.03);",
              div(
                class = "d-flex align-items-center justify-content-center gap-4 flex-wrap",
                tags$img(src = "scfigures/icar_logo.jpg", height = 90, width = 90, style = "border-radius: 10px; object-fit: contain;"),
                div(
                  class = "text-center text-md-start",
                  h4("Division of Agricultural Bioinformatics", style = "font-size: 1.35rem; color: #1E3A8A; margin-bottom: 3px; font-weight: 700;"),
                  h5("ICAR - Indian Agricultural Statistics Research Institute (IASRI)", style = "font-size: 1.10rem; color: #334155; margin-bottom: 3px; font-weight: 600;"),
                  p("New Delhi, India", style = "font-size: 0.95rem; color: #64748B; margin-bottom: 0; font-weight: 500;")
                )
              )
            ),
            
            # Author Team Grid: Row 1 (3 Authors)
            div(
              class = "row justify-content-center g-3 my-2",
              
              # 1. Kabilan Sakthivel
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/kabilan.JPG", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Kabilan Sakthivel", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Ph.D. Bioinformatics", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:kabilan151414@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "kabilan151414@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              ),
              
              # 2. Dr Dwijesh Chandra Mishra
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/mishra_sir.jpg", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Dr Dwijesh Chandra Mishra", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Principal Scientist", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:dwij.mishra@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "dwij.mishra@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              ),
              
              # 3. Dr Shashi Bhushan Lal
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/sb_lal_sir.JPG", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Dr Shashi Bhushan Lal", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Principal Scientist", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:sblall16@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "sblall16@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              )
            ),
            
            # Author Team Grid: Row 2 (3 Authors)
            div(
              class = "row justify-content-center g-3 my-2",
              
              # 4. Dr Sudhir Srivastava
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/sudhir_sir.JPG", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Dr Sudhir Srivastava", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Senior Scientist", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:sudhir0401bm@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "sudhir0401bm@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              ),
              
              # 5. Dr Krishna Kumar Chaturvedi
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/kkcsir.JPG", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Dr Krishna Kumar Chaturvedi", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Principal Scientist", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:kkcchaturvedi@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "kkcchaturvedi@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              ),
              
              # 6. Dr Sharanbasappa
              div(
                class = "col-md-4 col-sm-6 text-center",
                div(
                  class = "info-card", style = "background: #FFFFFF; border: 1px solid #E2E8F0; border-radius: 10px; padding: 16px 12px; height: 100%; box-shadow: 0 1px 3px rgba(0,0,0,0.03);",
                  tags$img(src = "scfigures/sharan_photo.jpg", height = 165, width = 135, style = "border-radius: 10px; border: 2px solid #E2E8F0; object-fit: cover;"),
                  p("Dr Sharanbasappa", style = "font-size: 0.96rem; margin-top: 10px; margin-bottom: 2px; font-weight: 700; color: #0F172A;"),
                  p("Scientist", style = "font-size: 0.82rem; margin-bottom: 6px; font-weight: 600; color: #0284C7;"),
                  p(
                    tags$a(
                      href = "mailto:smadival509@gmail.com",
                      style = "font-size: 0.78rem; color: #64748B; text-decoration: none;",
                      icon("envelope", class = "me-1 text-primary"),
                      "smadival509@gmail.com"
                    ),
                    style = "margin-bottom: 0;"
                  )
                )
              )
            ),
            
            # Feedback & GitHub Box
            div(
              style = "text-align: center; font-size: 0.88rem; color: #334155; margin: 26px auto 16px auto; max-width: 850px; padding: 14px 18px; background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 8px;",
              icon("github", class = "me-2 text-dark"),
              strong("For feedback, questions, or bug reports: "),
              "visit our GitHub repository at ",
              tags$a(href = "https://github.com/kabilanbio/scSimEval", target = "_blank", style = "color: #0284C7; text-decoration: none; font-weight: 600;", "github.com/kabilanbio/scSimEval")
            ),
            
            # Footer Note
            div(
              style = "font-size: 0.80rem; text-align: center; color: #64748B; background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px; margin-top: 10px;",
              span(icon("microscope", class = "me-1 text-primary"), tags$strong("Division of Agricultural Bioinformatics, ICAR-IASRI, New Delhi, India."))
            )
          )
        )
      )
    )
  ),
  nav_spacer(),
  nav_item(
    input_dark_mode(id = "dark_mode", mode = "light")
  )
)

# ==============================================================================
# SERVER DEFINITION
# ==============================================================================
server <- function(input, output, session) {
  
  # ============================================================================
  # LIVE RUNNING LOGS & DIAGNOSTIC ENGINE (v0.99.9)
  # ============================================================================
  # Session-level log storage container (plain environment for safe context-free reading/writing)
  log_storage <- new.env(parent = emptyenv())
  log_storage$entries <- list()
  log_storage$counter <- 0L
  
  # Reactive invalidation signal to notify Shiny outputs when new logs arrive
  sc_log_version <- reactiveVal(0L)

  sc_log <- function(level = c("INFO", "NOTE", "WARNING", "ERROR", "SUCCESS"), msg, component = "System", details = NULL) {
    level <- match.arg(level)
    cat(sprintf("[%s] [%-7s] [%-12s] %s\n", format(Sys.time(), "%H:%M:%S"), level, component, msg))
    log_storage$counter <- log_storage$counter + 1L
    new_entry <- list(
      id = log_storage$counter,
      time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      level = level,
      component = component,
      message = as.character(msg),
      details = if (!is.null(details) && nzchar(as.character(details))) as.character(details) else NULL
    )
    log_storage$entries[[log_storage$counter]] <- new_entry
    
    # Notify reactive consumers safely without error if outside reactive consumer
    tryCatch({
      isolate({
        sc_log_version(sc_log_version() + 1L)
      })
    }, error = function(e) NULL)
    
    invisible(new_entry)
  }

  # Startup Diagnostic Log Sequence
  sc_log("INFO", "scSimEval Interactive Benchmarking Studio v0.99.9 initialized.", "Startup")
  sc_log("NOTE", paste0("Runtime Environment: R ", getRversion(), " | Platform: ", R.version$platform, " | OS: ", Sys.info()[["sysname"]], " ", Sys.info()[["release"]]), "Environment")
  sc_log("INFO", paste0("System Concurrency: ", parallel::detectCores(), " logical CPU cores detected. Sparse matrix arithmetic (dgCMatrix) enabled."), "Hardware")
  sc_log("NOTE", "Evaluation Suite: 62 multiomics metrics registered across Categories 1 to 8.", "Metrics")
  
  # Check demo benchmark cache status
  demo_rds_path <- system.file("shiny/scSimEvalApp/demo_benchmark_data.rds", package = "scSimEval")
  if (nzchar(demo_rds_path) && file.exists(demo_rds_path)) {
    cache_sz <- round(file.info(demo_rds_path)$size / (1024 * 1024), 2)
    sc_log("INFO", paste0("Precomputed benchmark database found (", cache_sz, " MB). Fast demo loading enabled."), "Cache")
  } else {
    sc_log("NOTE", "Standalone demo cache file not pre-bundled; on-the-fly evaluation mode active.", "Cache")
  }
  sc_log("SUCCESS", "System ready. Navigate to 'Data Hub' to upload datasets or explore demo results.", "System")
  
  # Reactive state container
  # Reactive state container (Data loaded on-demand via options on sidebar)
  rv <- reactiveValues(
    benchmark_df = NULL,
    methods = NULL,
    toy_ref = NULL,
    toy_sim = NULL,
    sim_matrices = list(),
    cell_types = NULL,
    batch = NULL,
    source_name = "No Data Loaded",
    is_multidataset = FALSE,
    datasets = NULL,
    granular_benchmark_df = NULL,
    all_matrices = list(),
    ref_matrices = list(),
    dataset_summary_df = NULL
  )
  
  # Dedicated per-mode data cache ensuring complete isolation across upload modes
  mode_data <- reactiveValues(
    demo = list(
      benchmark_df = NULL, all_matrices = list(), ref_matrices = list(),
      sim_matrices = list(), dataset_summary_df = NULL, toy_ref = NULL,
      toy_sim = NULL, cell_types = NULL, batch = NULL, methods = NULL,
      source_name = "Demo Benchmark (Real Datasets)"
    ),
    scrna = list(
      benchmark_df = NULL, all_matrices = list(), ref_matrices = list(),
      sim_matrices = list(), dataset_summary_df = NULL, toy_ref = NULL,
      toy_sim = NULL, cell_types = NULL, batch = NULL, methods = NULL,
      source_name = "scRNA-seq Benchmarks"
    ),
    scatac = list(
      benchmark_df = NULL, all_matrices = list(), ref_matrices = list(),
      sim_matrices = list(), dataset_summary_df = NULL, toy_ref = NULL,
      toy_sim = NULL, cell_types = NULL, batch = NULL, methods = NULL,
      source_name = "scATAC-seq Benchmarks"
    ),
    unimodal = list(
      benchmark_df = NULL, all_matrices = list(), ref_matrices = list(),
      sim_matrices = list(), dataset_summary_df = NULL, toy_ref = NULL,
      toy_sim = NULL, cell_types = NULL, batch = NULL, methods = NULL,
      source_name = "User Evaluated Single-Cell Benchmarks"
    ),
    multiomics = list(
      benchmark_df = NULL, all_matrices = list(), ref_matrices = list(),
      sim_matrices = list(), dataset_summary_df = NULL, toy_ref = NULL,
      toy_sim = NULL, cell_types = NULL, batch = NULL, methods = NULL,
      source_name = "User Evaluated Multiomics Benchmarks"
    ),
    upload_bench = list(
      benchmark_df = NULL, all_matrices = list(), ref_matrices = list(),
      sim_matrices = list(), dataset_summary_df = NULL, toy_ref = NULL,
      toy_sim = NULL, cell_types = NULL, batch = NULL, methods = NULL,
      source_name = "Uploaded Benchmark Results"
    )
  )

  # Synchronize active reactive values (rv) with active data mode
  sync_active_mode_data <- function(target_mode = NULL) {
    target_mode <- target_mode %||% (input$opt_data_mode %||% "demo")
    cur <- mode_data[[target_mode]]
    if (!is.null(cur) && (!is.null(cur$benchmark_df) || length(cur$all_matrices) > 0 || !is.null(cur$dataset_summary_df))) {
      rv$benchmark_df       <- cur$benchmark_df
      rv$methods            <- cur$methods
      rv$all_matrices       <- cur$all_matrices
      rv$ref_matrices       <- cur$ref_matrices
      rv$sim_matrices       <- cur$sim_matrices
      rv$dataset_summary_df <- cur$dataset_summary_df
      rv$toy_ref            <- cur$toy_ref
      rv$toy_sim            <- cur$toy_sim
      rv$cell_types         <- cur$cell_types
      rv$batch              <- cur$batch
      rv$cell_types_map     <- cur$cell_types_map
      rv$batches_map        <- cur$batches_map
      rv$source_name        <- cur$source_name
      if (!is.null(cur$methods) && length(cur$methods) > 0) {
        updateCheckboxGroupInput(session, "sel_bubble_methods", choices = cur$methods, selected = cur$methods)
      }
    } else {
      rv$benchmark_df       <- NULL
      rv$methods            <- NULL
      rv$all_matrices       <- list()
      rv$ref_matrices       <- list()
      rv$sim_matrices       <- list()
      rv$dataset_summary_df <- NULL
      rv$toy_ref            <- NULL
      rv$toy_sim            <- NULL
      rv$cell_types         <- NULL
      rv$batch              <- NULL
      rv$cell_types_map     <- list()
      rv$batches_map        <- list()
      rv$source_name        <- paste0("No data loaded for ", target_mode)
    }
  }

  observeEvent(input$opt_data_mode, {
    sync_active_mode_data(input$opt_data_mode)
  }, ignoreInit = TRUE)
  
  # Reactive tracking for user-confirmed plotting (prevents unwanted auto-plotting)
  rv_plot_triggers <- reactiveValues(
    bubble = 0,
    sum    = 0,
    dist   = 0,
    scale  = 0,
    box    = 0,
    heat   = 0,
    pca    = 0,
    mds    = 0,
    emb    = 0
  )
  
  # Helper for clean confirmation placeholder
  render_plot_placeholder <- function(title, desc, button_id, btn_label = "Confirm & Generate Plot", icon_name = "chart-simple") {
    div(
      class = "plot-placeholder-card",
      icon(icon_name, class = "text-primary mb-3", style = "font-size: 2.6rem; opacity: 0.9;"),
      h5(title, class = "placeholder-title"),
      p(desc, class = "placeholder-desc"),
      actionButton(button_id, btn_label, class = "btn btn-primary px-4 py-2", icon = icon("play"), style = "font-weight: 600; box-shadow: 0 2px 6px rgba(30,58,138,0.25);")
    )
  }
  
  # Plot confirmation observers for all 9 plots
  observeEvent(input$btn_plot_bubble, { rv_plot_triggers$bubble <- rv_plot_triggers$bubble + 1 })
  observeEvent(input$btn_plot_bubble_inline, { rv_plot_triggers$bubble <- rv_plot_triggers$bubble + 1 })
  observeEvent(input$btn_plot_sum, { rv_plot_triggers$sum <- rv_plot_triggers$sum + 1 })
  observeEvent(input$btn_plot_sum_inline, { rv_plot_triggers$sum <- rv_plot_triggers$sum + 1 })
  observeEvent(input$btn_plot_dist, { rv_plot_triggers$dist <- rv_plot_triggers$dist + 1 })
  observeEvent(input$btn_plot_dist_inline, { rv_plot_triggers$dist <- rv_plot_triggers$dist + 1 })
  observeEvent(input$btn_plot_scale, { rv_plot_triggers$scale <- rv_plot_triggers$scale + 1 })
  observeEvent(input$btn_plot_scale_inline, { rv_plot_triggers$scale <- rv_plot_triggers$scale + 1 })
  observeEvent(input$btn_plot_box, { rv_plot_triggers$box <- rv_plot_triggers$box + 1 })
  observeEvent(input$btn_plot_box_inline, { rv_plot_triggers$box <- rv_plot_triggers$box + 1 })
  observeEvent(input$btn_plot_heat, { rv_plot_triggers$heat <- rv_plot_triggers$heat + 1 })
  observeEvent(input$btn_plot_heat_inline, { rv_plot_triggers$heat <- rv_plot_triggers$heat + 1 })
  observeEvent(input$btn_plot_pca, { rv_plot_triggers$pca <- rv_plot_triggers$pca + 1 })
  observeEvent(input$btn_plot_pca_inline, { rv_plot_triggers$pca <- rv_plot_triggers$pca + 1 })
  observeEvent(input$btn_plot_mds, { rv_plot_triggers$mds <- rv_plot_triggers$mds + 1 })
  observeEvent(input$btn_plot_mds_inline, { rv_plot_triggers$mds <- rv_plot_triggers$mds + 1 })
  observeEvent(input$btn_plot_emb, { rv_plot_triggers$emb <- rv_plot_triggers$emb + 1 })
  observeEvent(input$btn_plot_emb_inline, { rv_plot_triggers$emb <- rv_plot_triggers$emb + 1 })
  
  # Navigation triggers
  observeEvent(input$btn_go_data, { nav_select("nav_active", "Data Hub") })
  observeEvent(input$btn_go_data_inline, { nav_select("nav_active", "Data Hub") })
  observeEvent(input$btn_go_bubble, { nav_select("nav_active", "Comparative Bubble Matrix") })
  observeEvent(input$btn_go_viz, { nav_select("nav_active", "Visualizations") })
  observeEvent(input$btn_go_download, { nav_select("nav_active", "Benchmark Results & Export") })
  observeEvent(input$btn_go_help, { nav_select("nav_active", "Help & Getting Started") })
  observeEvent(input$btn_go_contact, { nav_select("nav_active", "Contact") })
  observeEvent(input$btn_empty_help, { nav_select("nav_active", "Help & Getting Started") })
  # Modal format guide for cell types and batches
  show_metadata_format_modal <- function() {
    showModal(modalDialog(
      title = div(icon("table-list", class = "me-2 text-primary"), tags$strong("Metadata Input Formats (Cell Types & Batches)")),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      div(
        p("Cell-type and batch annotations are optional inputs used to compute clustering metrics (ARI, NMI, Silhouette), marker gene DEGs, and batch mixing confounder alignment:"),
        # Metric measurement capacity callout
        div(
          style = "background: #F0FDF4; border: 1px solid #BBF7D0; border-radius: 6px; padding: 14px 16px; margin-bottom: 16px;",
          h6(tags$b(icon("chart-line", class = "me-1 text-success"), " Quantitative Metric Measurement Capacity by Upload Configuration"), style = "color: #166534; font-size: 0.94rem; margin-bottom: 6px;"),
          p("Depending on whether metadata files are provided, scSimEval scales the quantitative evaluation suite dynamically across categories:", style = "font-size: 0.82rem; color: #15803D; margin-bottom: 10px;"),
          tags$div(
            class = "table-responsive",
            tags$table(
              class = "table table-sm table-bordered mb-1 bg-white align-middle",
              style = "font-size: 0.79rem;",
              tags$thead(
                class = "table-light",
                tags$tr(
                  tags$th(style = "width: 18%;", "Experiment Mode"),
                  tags$th(style = "width: 25%;", "Uploaded Inputs"),
                  tags$th(style = "width: 17%; text-align: center;", "Metrics Measured"),
                  tags$th(style = "width: 40%;", "Active Evaluation Categories & Scope")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td(tags$b("1. scRNA-seq")),
                  tags$td("Count Matrix alone (No metadata)"),
                  tags$td(style = "text-align: center;", span(class = "badge", style = "background-color: #475569 !important; color: #FFFFFF !important; font-weight: 700; font-size: 0.82rem; padding: 5px 10px; border-radius: 6px;", "24 / 62 Metrics")),
                  tags$td("Cat I (Distributional Properties, 14 moments & 2D tests), Cat II (Correlations & Zero-Inflation, 6 metrics), Cat VI (Trajectory dynamics, 2 metrics), Cat VIII (Computational Scalability, 2 metrics)")
                ),
                tags$tr(
                  tags$td(tags$b("1. scRNA-seq")),
                  tags$td("Count Matrix + Cell Types"),
                  tags$td(style = "text-align: center;", span(class = "badge", style = "background-color: #0284C7 !important; color: #FFFFFF !important; font-weight: 700; font-size: 0.82rem; padding: 5px 10px; border-radius: 6px;", "44 / 62 Metrics")),
                  tags$td("Adds +20 metrics: Unlocks Cat III (Cellular Structure & Clustering Concordance: ARI, NMI, AMI, Silhouette, Dunn, DB, Calinski-Harabasz, Neighborhood Purity, 10 metrics) and Cat V (Biological Signal Fidelity: SimBench DE, DV, DD, DP, top DEG Jaccard, delta variance, 10 metrics)")
                ),
                tags$tr(
                  tags$td(tags$b("1. scRNA-seq")),
                  tags$td(tags$b("Full Metadata (Matrix + Cell Types + Batch)")),
                  tags$td(style = "text-align: center;", span(class = "badge", style = "background-color: #16A34A !important; color: #FFFFFF !important; font-weight: 700; font-size: 0.82rem; padding: 5px 10px; border-radius: 6px;", "51 / 62 Metrics")),
                  tags$td("Adds +7 metrics: Unlocks Cat IV (Batch Effects & Confounder Mixing: Batch Silhouette, Shannon entropy, PCR R², CMS, ISI, Seurat mixing, LDF difference). Maximum possible for unimodal RNA.")
                ),
                tags$tr(
                  tags$td(tags$b("2. scATAC-seq")),
                  tags$td("Count Matrix alone (No metadata)"),
                  tags$td(style = "text-align: center;", span(class = "badge", style = "background-color: #475569 !important; color: #FFFFFF !important; font-weight: 700; font-size: 0.82rem; padding: 5px 10px; border-radius: 6px;", "22 / 62 Metrics")),
                  tags$td("Cat I (Peak accessibility distribution moments, 14 metrics), Cat II (Peak-peak correlations & sparsity curves, 6 metrics), Cat VIII (Computational Scalability, 2 metrics)")
                ),
                tags$tr(
                  tags$td(tags$b("2. scATAC-seq")),
                  tags$td("Count Matrix + Cell Types"),
                  tags$td(style = "text-align: center;", span(class = "badge", style = "background-color: #0284C7 !important; color: #FFFFFF !important; font-weight: 700; font-size: 0.82rem; padding: 5px 10px; border-radius: 6px;", "36 / 62 Metrics")),
                  tags$td("Adds +14 metrics: Unlocks Cat III (Epigenomic Cellular Clustering Concordance, 8 metrics) and Cat V (Differential Accessibility signal fidelity, 6 metrics)")
                ),
                tags$tr(
                  tags$td(tags$b("2. scATAC-seq")),
                  tags$td(tags$b("Full Metadata (Matrix + Cell Types + Batch)")),
                  tags$td(style = "text-align: center;", span(class = "badge", style = "background-color: #16A34A !important; color: #FFFFFF !important; font-weight: 700; font-size: 0.82rem; padding: 5px 10px; border-radius: 6px;", "43 / 62 Metrics")),
                  tags$td("Adds +7 metrics: Unlocks Cat IV (Epigenomic Batch Effects & Confounder Mixing: Batch Silhouette, Shannon entropy, PCR R², CMS, ISI). Maximum possible for unimodal ATAC.")
                ),
                tags$tr(
                  tags$td(tags$b("3. Paired Multiomics")),
                  tags$td("RNA + ATAC Matrices alone (No metadata)"),
                  tags$td(style = "text-align: center;", span(class = "badge", style = "background-color: #4F46E5 !important; color: #FFFFFF !important; font-weight: 700; font-size: 0.82rem; padding: 5px 10px; border-radius: 6px;", "36 / 62 Metrics")),
                  tags$td("Measures both unimodal layers simultaneously + Cat VII Cross-Modal Coupling (FOSCTTM cell alignment, Match@1, Match@5, cross-modal generation, RV coefficient, peak-to-gene linkage, module co-expression, 10 metrics)")
                ),
                tags$tr(
                  tags$td(tags$b("3. Paired Multiomics")),
                  tags$td("RNA + ATAC + Cell Types"),
                  tags$td(style = "text-align: center;", span(class = "badge", style = "background-color: #0284C7 !important; color: #FFFFFF !important; font-weight: 700; font-size: 0.82rem; padding: 5px 10px; border-radius: 6px;", "55 / 62 Metrics")),
                  tags$td("Unlocks Cat III (Dual cellular clustering concordance), Cat V (Dual DEG/DA signal fidelity), and Cat VII (Cell-type guided cross-modal alignment accuracy & Macro-F1)")
                ),
                tags$tr(
                  tags$td(tags$b("3. Paired Multiomics")),
                  tags$td(tags$b("Full Metadata (RNA + ATAC + Cell Types + Batch)")),
                  tags$td(style = "text-align: center;", span(class = "badge", style = "background-color: #15803D !important; color: #FFFFFF !important; font-weight: 700; font-size: 0.82rem; padding: 5px 10px; border-radius: 6px;", "All 62 / 62 Metrics")),
                  tags$td("100% complete coverage across all 8 biological and computational dimensions (Categories I–VIII)")
                ),
                tags$tr(
                  tags$td(tags$b("4. Unpaired Multiomics")),
                  tags$td("RNA + ATAC Matrices alone (Independent cells)"),
                  tags$td(style = "text-align: center;", span(class = "badge", style = "background-color: #475569 !important; color: #FFFFFF !important; font-weight: 700; font-size: 0.82rem; padding: 5px 10px; border-radius: 6px;", "30 / 62 Metrics")),
                  tags$td("Unimodal distributions for both layers + population-level co-expression module fidelity, Network Jaccard, and cross-modality correlation")
                ),
                tags$tr(
                  tags$td(tags$b("4. Unpaired Multiomics")),
                  tags$td(tags$b("RNA + ATAC + Cell Types (+ Batch)")),
                  tags$td(style = "text-align: center;", span(class = "badge", style = "background-color: #16A34A !important; color: #FFFFFF !important; font-weight: 700; font-size: 0.82rem; padding: 5px 10px; border-radius: 6px;", "52 / 62 Metrics")),
                  tags$td("Unlocks cluster-matched cross-modal label transfer accuracy, cell-type regulatory concordance, and technical donor batch alignment")
                )
              )
            )
          )
        ),

        fluidRow(
          column(6,
            div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 14px; margin-bottom: 12px;",
              h6(tags$b(icon("layer-group", class = "me-1 text-primary"), " 1. Cell Types Format"), style = "color: #1E3A8A; margin-bottom: 6px; font-size: 0.90rem;"),
              p("Upload a 2-column CSV/TSV table with ", tags$code("Cell_ID"), " and ", tags$code("Cell_Type"), " (or a 1-column vector in the same cell order as the count matrix):", style = "font-size: 0.80rem; color: #475569; margin-bottom: 8px;"),
              tags$pre("Cell_ID,Cell_Type\nCell_01,TypeA\nCell_02,TypeA\nCell_03,TypeB\nCell_04,TypeB\n...", style = "font-size: 0.76rem; background: #FFFFFF; border: 1px solid #CBD5E1; padding: 8px; border-radius: 4px; margin-bottom: 10px;"),
              downloadButton("dl_modal_celltypes_csv", "Download Example Cell Types (.csv)", class = "btn btn-sm btn-outline-primary w-100", icon = icon("download"))
            )
          ),
          column(6,
            div(style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 14px; margin-bottom: 12px;",
              h6(tags$b(icon("boxes-stacked", class = "me-1 text-teal"), " 2. Batch Annotations Format"), style = "color: #0D9488; margin-bottom: 6px; font-size: 0.90rem;"),
              p("Upload a 2-column CSV/TSV table with ", tags$code("Cell_ID"), " and ", tags$code("Batch"), " (or a 1-column vector in the same cell order as the count matrix):", style = "font-size: 0.80rem; color: #475569; margin-bottom: 8px;"),
              tags$pre("Cell_ID,Batch\nCell_01,Batch1\nCell_02,Batch2\nCell_03,Batch1\nCell_04,Batch2\n...", style = "font-size: 0.76rem; background: #FFFFFF; border: 1px solid #CBD5E1; padding: 8px; border-radius: 4px; margin-bottom: 10px;"),
              downloadButton("dl_modal_batch_csv", "Download Example Batch (.csv)", class = "btn btn-sm btn-outline-success w-100", icon = icon("download"))
            )
          )
        )
      )
    ))
  }
  observeEvent(input$btn_show_metadata_format, { show_metadata_format_modal() })
  observeEvent(input$btn_show_metadata_format_multi, { show_metadata_format_modal() })
  observeEvent(input$btn_show_metadata_format_hint, { show_metadata_format_modal() })
  observeEvent(input$btn_show_metadata_format_multi_hint, { show_metadata_format_modal() })

  # Helper function to generate sample metadata data.frame
  get_sample_celltypes_df <- function() {
    data(example_scrna, package = "scSimEval")
    data.frame(
      Cell_ID = names(example_scrna$cell_types),
      Cell_Type = as.character(example_scrna$cell_types),
      stringsAsFactors = FALSE
    )
  }
  get_sample_batch_df <- function() {
    data(example_scrna, package = "scSimEval")
    data.frame(
      Cell_ID = names(example_scrna$batch_info),
      Batch = as.character(example_scrna$batch_info),
      stringsAsFactors = FALSE
    )
  }

  # Download Handlers for sample metadata files
  output$dl_sample_celltypes_csv <- downloadHandler(
    filename = function() { "example_cell_types.csv" },
    content = function(file) { utils::write.csv(get_sample_celltypes_df(), file, row.names = FALSE) }
  )
  output$dl_sample_celltypes_csv_multi <- downloadHandler(
    filename = function() { "example_cell_types.csv" },
    content = function(file) { utils::write.csv(get_sample_celltypes_df(), file, row.names = FALSE) }
  )
  output$dl_modal_celltypes_csv <- downloadHandler(
    filename = function() { "example_cell_types.csv" },
    content = function(file) { utils::write.csv(get_sample_celltypes_df(), file, row.names = FALSE) }
  )
  output$dl_help_sample_celltypes_csv <- downloadHandler(
    filename = function() { "example_cell_types.csv" },
    content = function(file) { utils::write.csv(get_sample_celltypes_df(), file, row.names = FALSE) }
  )
  output$dl_modal_batch_csv <- downloadHandler(
    filename = function() { "example_batch_annotations.csv" },
    content = function(file) { utils::write.csv(get_sample_batch_df(), file, row.names = FALSE) }
  )
  output$dl_help_sample_batch_csv <- downloadHandler(
    filename = function() { "example_batch_annotations.csv" },
    content = function(file) { utils::write.csv(get_sample_batch_df(), file, row.names = FALSE) }
  )
  
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 1 - Load Demo Benchmark
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # Data Hub: Mode 1 - Load Example Dataset Analysis
  # ----------------------------------------------------------------------------
  # Realistic baseline scalability helper across simulators and datasets
  # Realistic baseline scalability helper across simulators and datasets
  get_default_sim_scalability <- function(s, d, mode = 'unimodal') {
    base_time_uni <- c(12.4, 18.6, 8.2, 15.1, 11.5, 14.8)
    base_mem_uni  <- c(248, 312, 185, 290, 220, 275)
    d_scale_time <- if (d == 1) 1.0 else (1.0 + (d - 1) * 0.22)
    d_scale_mem  <- if (d == 1) 1.0 else (1.0 + (d - 1) * 0.13)
    mult <- if (identical(mode, 'multiomics')) 1.4 else 1.0
    idx <- ((s - 1) %% 6) + 1
    t_val <- round(base_time_uni[idx] * d_scale_time * mult, 1)
    m_val <- round(base_mem_uni[idx] * d_scale_mem * mult, 0)
    list(time = t_val, mem = m_val)
  }
  
  # Helper to locate package data files across installation & development roots
  get_pkg_data_file <- function(filename) {
    candidates <- c(
      system.file("shiny", "scSimEvalApp", "data", filename, package = "scSimEval"),
      system.file("data", filename, package = "scSimEval"),
      file.path("data", filename),
      file.path("inst", "shiny", "scSimEvalApp", "data", filename),
      file.path("..", "data", filename),
      file.path("..", "..", "inst", "shiny", "scSimEvalApp", "data", filename)
    )
    for (f in candidates) {
      if (nzchar(f) && file.exists(f)) return(normalizePath(f, winslash = "/"))
    }
    return(NULL)
  }

  # Unified Example Dataset & Benchmark Loader (Real Rice Bud & Soybean Root Datasets)
  load_example_benchmark <- function(type = c("all", "scrna_rice", "scrna_soy", "scatac_rice", "scatac_soy", "multi_paired", "multi_unpaired")) {
    type <- match.arg(type)
    sc_log("INFO", paste0("Requesting demo benchmark data for configuration: '", type, "'..."), "DataHub")
    
    bm_file <- get_pkg_data_file("demo_benchmark_data.rds")
    rb_file <- get_pkg_data_file("example_rice_bud.rds")
    sb_file <- get_pkg_data_file("example_soybean_root.rds")
    
    if (is.null(bm_file) || !file.exists(bm_file)) {
      sc_log("ERROR", "Benchmark database cache file demo_benchmark_data.rds not found on disk.", "DataHub")
      showNotification("Benchmark database cache file not found on disk.", type = "error")
      return(FALSE)
    }
    
    bm <- readRDS(bm_file)
    rb <- if (!is.null(rb_file) && file.exists(rb_file)) readRDS(rb_file) else NULL
    sb <- if (!is.null(sb_file) && file.exists(sb_file)) readRDS(sb_file) else NULL
    
    get_rb_ct <- function() {
      if (is.null(rb) || is.null(rb$reference$cell_types)) return(NULL)
      df <- rb$reference$cell_types
      ct <- df$CellType %||% (if (ncol(df) >= 2) df[[2]] else df[[1]])
      as.character(ct)
    }
    get_rb_bt <- function() {
      if (is.null(rb) || is.null(rb$reference$batches)) return(NULL)
      df <- rb$reference$batches
      bt <- df$Batch %||% (if (ncol(df) >= 2) df[[2]] else df[[1]])
      as.character(bt)
    }
    get_sb_ct <- function(mod = "rna") {
      if (is.null(sb)) return(NULL)
      df <- if (tolower(mod) == "atac") sb$reference$atac_cell_types else sb$reference$rna_cell_types
      if (is.null(df)) return(NULL)
      ct <- df$CellType %||% (if (ncol(df) >= 2) df[[2]] else df[[1]])
      as.character(ct)
    }
    get_sb_bt <- function(mod = "rna") {
      if (is.null(sb)) return(NULL)
      df <- if (tolower(mod) == "atac") sb$reference$atac_batches else sb$reference$rna_batches
      if (is.null(df)) return(NULL)
      bt <- df$Batch %||% (if (ncol(df) >= 2) df[[2]] else df[[1]])
      as.character(bt)
    }
    
    target_mode <- "multiomics"
    ref_mats <- list()
    sim_mats <- list()
    all_mats <- list()
    ct <- NULL
    bt <- NULL
    bench_df <- NULL
    summary_df <- NULL
    source_name <- ""
    
    if (type == "scrna_rice") {
      target_mode <- "scrna"
      req(rb)
      ref_mats <- list("Rice Bud (scRNA-seq)" = rb$reference$rna)
      sim_mats <- lapply(rb$simulated, function(x) x$rna)
      all_mats <- c(list("[Reference] Rice Bud (scRNA-seq)" = rb$reference$rna),
                    stats::setNames(sim_mats, paste0("[Simulated] ", names(sim_mats), " (Rice Bud scRNA-seq)")))
      bench_df <- bm$benchmark_df[bm$benchmark_df$Dataset == "Rice Bud (Paired)" & bm$benchmark_df$Modality %in% c("RNA", "Computation"), ]
      bench_df$Dataset <- "Rice Bud (scRNA-seq)"
      summary_df <- bm$dataset_summary_df[bm$dataset_summary_df$Dataset == "Rice Bud (Paired)" & bm$dataset_summary_df$Modality == "scRNA-seq", ]
      summary_df$Dataset <- "Rice Bud (scRNA-seq)"
      ct <- get_rb_ct()
      bt <- get_rb_bt()
      source_name <- "Example Data: Rice Bud scRNA-seq (100 cells, 200 genes, 8 cell types, 6 simulators)"
      
    } else if (type == "scrna_soy") {
      target_mode <- "scrna"
      req(sb)
      ref_mats <- list("Soybean Root (scRNA-seq)" = sb$reference$rna)
      sim_mats <- lapply(sb$simulated, function(x) x$rna)
      all_mats <- c(list("[Reference] Soybean Root (scRNA-seq)" = sb$reference$rna),
                    stats::setNames(sim_mats, paste0("[Simulated] ", names(sim_mats), " (Soybean Root scRNA-seq)")))
      bench_df <- bm$benchmark_df[bm$benchmark_df$Dataset == "Soybean Root (Unpaired)" & bm$benchmark_df$Modality %in% c("RNA", "Computation"), ]
      bench_df$Dataset <- "Soybean Root (scRNA-seq)"
      summary_df <- bm$dataset_summary_df[bm$dataset_summary_df$Dataset == "Soybean Root (Unpaired)" & bm$dataset_summary_df$Modality == "scRNA-seq", ]
      summary_df$Dataset <- "Soybean Root (scRNA-seq)"
      ct <- get_sb_ct("rna")
      bt <- get_sb_bt("rna")
      source_name <- "Example Data: Soybean Root scRNA-seq (100 cells, 200 genes, 8 cell types, 2 batches, 3 simulators)"
      
    } else if (type == "scatac_rice") {
      target_mode <- "scatac"
      req(rb)
      ref_mats <- list("Rice Bud (scATAC-seq)" = rb$reference$atac)
      sim_mats <- lapply(rb$simulated, function(x) x$atac)
      all_mats <- c(list("[Reference] Rice Bud (scATAC-seq)" = rb$reference$atac),
                    stats::setNames(sim_mats, paste0("[Simulated] ", names(sim_mats), " (Rice Bud scATAC-seq)")))
      bench_df <- bm$benchmark_df[bm$benchmark_df$Dataset == "Rice Bud (Paired)" & bm$benchmark_df$Modality %in% c("ATAC", "Computation"), ]
      bench_df$Dataset <- "Rice Bud (scATAC-seq)"
      summary_df <- bm$dataset_summary_df[bm$dataset_summary_df$Dataset == "Rice Bud (Paired)" & bm$dataset_summary_df$Modality == "scATAC-seq", ]
      summary_df$Dataset <- "Rice Bud (scATAC-seq)"
      ct <- get_rb_ct()
      bt <- get_rb_bt()
      source_name <- "Example Data: Rice Bud scATAC-seq (100 cells, 200 peaks, 8 cell types, 6 simulators)"
      
    } else if (type == "scatac_soy") {
      target_mode <- "scatac"
      req(sb)
      ref_mats <- list("Soybean Root (scATAC-seq)" = sb$reference$atac)
      sim_mats <- lapply(sb$simulated, function(x) x$atac)
      all_mats <- c(list("[Reference] Soybean Root (scATAC-seq)" = sb$reference$atac),
                    stats::setNames(sim_mats, paste0("[Simulated] ", names(sim_mats), " (Soybean Root scATAC-seq)")))
      bench_df <- bm$benchmark_df[bm$benchmark_df$Dataset == "Soybean Root (Unpaired)" & bm$benchmark_df$Modality %in% c("ATAC", "Computation"), ]
      bench_df$Dataset <- "Soybean Root (scATAC-seq)"
      summary_df <- bm$dataset_summary_df[bm$dataset_summary_df$Dataset == "Soybean Root (Unpaired)" & bm$dataset_summary_df$Modality == "scATAC-seq", ]
      summary_df$Dataset <- "Soybean Root (scATAC-seq)"
      ct <- get_sb_ct("atac")
      bt <- get_sb_bt("atac")
      source_name <- "Example Data: Soybean Root scATAC-seq (100 cells, 200 peaks, 10 cell types, 5 batches, 3 simulators)"
      
    } else if (type == "multi_paired") {
      target_mode <- "multiomics"
      req(rb)
      ref_mats <- list("Rice Bud_RNA" = rb$reference$rna, "Rice Bud_ATAC" = rb$reference$atac)
      all_mats <- list(
        "[Reference] Rice Bud RNA" = rb$reference$rna,
        "[Reference] Rice Bud ATAC" = rb$reference$atac
      )
      for (m in names(rb$simulated)) {
        sim_mats[[paste0(m, "_RNA_Rice Bud")]] <- rb$simulated[[m]]$rna
        sim_mats[[paste0(m, "_ATAC_Rice Bud")]] <- rb$simulated[[m]]$atac
        all_mats[[paste0("[Simulated] ", m, " RNA (Rice Bud)")]] <- rb$simulated[[m]]$rna
        all_mats[[paste0("[Simulated] ", m, " ATAC (Rice Bud)")]] <- rb$simulated[[m]]$atac
      }
      bench_df <- bm$benchmark_df[bm$benchmark_df$Dataset == "Rice Bud (Paired)", ]
      summary_df <- bm$dataset_summary_df[bm$dataset_summary_df$Dataset == "Rice Bud (Paired)", ]
      ct <- get_rb_ct()
      bt <- get_rb_bt()
      source_name <- "Example Data: Rice Bud Paired Multiomics (100 cells, 200 genes, 200 peaks, 8 cell types, 6 simulators)"
      
    } else if (type == "multi_unpaired") {
      target_mode <- "multiomics"
      req(sb)
      ref_mats <- list("Soybean Root_RNA" = sb$reference$rna, "Soybean Root_ATAC" = sb$reference$atac)
      all_mats <- list(
        "[Reference] Soybean Root RNA" = sb$reference$rna,
        "[Reference] Soybean Root ATAC" = sb$reference$atac
      )
      for (m in names(sb$simulated)) {
        sim_mats[[paste0(m, "_RNA_Soybean Root")]] <- sb$simulated[[m]]$rna
        sim_mats[[paste0(m, "_ATAC_Soybean Root")]] <- sb$simulated[[m]]$atac
        all_mats[[paste0("[Simulated] ", m, " RNA (Soybean Root)")]] <- sb$simulated[[m]]$rna
        all_mats[[paste0("[Simulated] ", m, " ATAC (Soybean Root)")]] <- sb$simulated[[m]]$atac
      }
      bench_df <- bm$benchmark_df[bm$benchmark_df$Dataset == "Soybean Root (Unpaired)", ]
      summary_df <- bm$dataset_summary_df[bm$dataset_summary_df$Dataset == "Soybean Root (Unpaired)", ]
      ct <- get_sb_ct("rna")
      bt <- get_sb_bt("rna")
      source_name <- "Example Data: Soybean Root Unpaired Multiomics (100 RNA cells, 100 ATAC cells, 200 genes, 200 peaks, 3 simulators)"
      
    } else {
      # "all" - Full benchmark of both datasets
      target_mode <- input$opt_data_mode %||% "multiomics"
      if (!target_mode %in% c("scrna", "scatac", "multiomics")) target_mode <- "multiomics"
      bench_df <- bm$benchmark_df
      summary_df <- bm$dataset_summary_df
      if (!is.null(rb)) {
        ref_mats[["Rice Bud_RNA"]] <- rb$reference$rna
        ref_mats[["Rice Bud_ATAC"]] <- rb$reference$atac
        all_mats[["[Reference] Rice Bud RNA"]] <- rb$reference$rna
        all_mats[["[Reference] Rice Bud ATAC"]] <- rb$reference$atac
        for (m in names(rb$simulated)) {
          sim_mats[[paste0(m, "_RNA_Rice Bud")]] <- rb$simulated[[m]]$rna
          sim_mats[[paste0(m, "_ATAC_Rice Bud")]] <- rb$simulated[[m]]$atac
          all_mats[[paste0("[Simulated] ", m, " RNA (Rice Bud)")]] <- rb$simulated[[m]]$rna
          all_mats[[paste0("[Simulated] ", m, " ATAC (Rice Bud)")]] <- rb$simulated[[m]]$atac
        }
        ct <- get_rb_ct()
        bt <- get_rb_bt()
      }
      if (!is.null(sb)) {
        ref_mats[["Soybean Root_RNA"]] <- sb$reference$rna
        ref_mats[["Soybean Root_ATAC"]] <- sb$reference$atac
        all_mats[["[Reference] Soybean Root RNA"]] <- sb$reference$rna
        all_mats[["[Reference] Soybean Root ATAC"]] <- sb$reference$atac
        for (m in names(sb$simulated)) {
          sim_mats[[paste0(m, "_RNA_Soybean Root")]] <- sb$simulated[[m]]$rna
          sim_mats[[paste0(m, "_ATAC_Soybean Root")]] <- sb$simulated[[m]]$atac
          all_mats[[paste0("[Simulated] ", m, " RNA (Soybean Root)")]] <- sb$simulated[[m]]$rna
          all_mats[[paste0("[Simulated] ", m, " ATAC (Soybean Root)")]] <- sb$simulated[[m]]$atac
        }
      }
      source_name <- "Benchmark Database: Rice Bud (6 Simulators) & Soybean Root (3 Simulators)"
    }
    
    methods <- unique(bench_df$Method)
    
    ct_map <- list()
    bt_map <- list()
    if (!is.null(rb)) {
      rb_ct <- get_rb_ct()
      rb_bt <- get_rb_bt()
      ct_map[["Rice Bud_RNA"]] <- rb_ct
      ct_map[["Rice Bud_ATAC"]] <- rb_ct
      ct_map[["Rice Bud"]] <- rb_ct
      ct_map[["Rice Bud (scRNA-seq)_RNA"]] <- rb_ct
      ct_map[["Rice Bud (scATAC-seq)_ATAC"]] <- rb_ct
      bt_map[["Rice Bud_RNA"]] <- rb_bt
      bt_map[["Rice Bud_ATAC"]] <- rb_bt
      bt_map[["Rice Bud"]] <- rb_bt
    }
    if (!is.null(sb)) {
      sb_ct_r <- get_sb_ct("rna")
      sb_ct_a <- get_sb_ct("atac")
      sb_bt_r <- get_sb_bt("rna")
      sb_bt_a <- get_sb_bt("atac")
      ct_map[["Soybean Root_RNA"]] <- sb_ct_r
      ct_map[["Soybean Root_ATAC"]] <- sb_ct_a
      ct_map[["Soybean Root"]] <- sb_ct_r
      ct_map[["Soybean Root (scRNA-seq)_RNA"]] <- sb_ct_r
      ct_map[["Soybean Root (scATAC-seq)_ATAC"]] <- sb_ct_a
      bt_map[["Soybean Root_RNA"]] <- sb_bt_r
      bt_map[["Soybean Root_ATAC"]] <- sb_bt_a
      bt_map[["Soybean Root"]] <- sb_bt_r
    }
    
    mode_data[[target_mode]] <- list(
      benchmark_df       = bench_df,
      methods            = methods,
      toy_ref            = if (length(ref_mats) > 0) ref_mats[[1]] else NULL,
      toy_sim            = if (length(sim_mats) > 0) sim_mats[[1]] else NULL,
      sim_matrices       = sim_mats,
      ref_matrices       = ref_mats,
      cell_types         = ct,
      batch              = bt,
      cell_types_map     = ct_map,
      batches_map        = bt_map,
      source_name        = source_name,
      dataset_summary_df = summary_df,
      all_matrices       = all_mats
    )
    mode_data$demo <- mode_data[[target_mode]]
    
    sync_active_mode_data(target_mode)
    rv_plot_triggers$bubble <- 0; rv_plot_triggers$sum <- 0; rv_plot_triggers$dist <- 0; rv_plot_triggers$scale <- 0; rv_plot_triggers$box <- 0; rv_plot_triggers$heat <- 0; rv_plot_triggers$pca <- 0; rv_plot_triggers$mds <- 0; rv_plot_triggers$emb <- 0
    sc_log("SUCCESS", paste0("Demo benchmark data loaded: '", source_name, "' (", nrow(bench_df), " evaluation rows across ", length(methods), " simulators: ", paste(methods, collapse = ", "), ")."), "DataHub")
    showNotification(paste0(source_name, " loaded successfully!"), type = "message", duration = 6)
    nav_select("datahub_subtabs", "Dataset properties")
    return(TRUE)
  }
  
  execute_load_demo_benchmark <- function() {
    load_example_benchmark("all")
  }
  
  observeEvent(input$btn_load_example_scrna, {
    ds <- input$opt_example_scrna_dataset %||% "rice_bud"
    load_example_benchmark(if (ds == "soybean_root") "scrna_soy" else "scrna_rice")
  })
  observeEvent(input$btn_load_example_scatac, {
    ds <- input$opt_example_scatac_dataset %||% "rice_bud"
    load_example_benchmark(if (ds == "soybean_root") "scatac_soy" else "scatac_rice")
  })
  observeEvent(input$btn_load_example_paired_multi, {
    load_example_benchmark("multi_paired")
  })
  observeEvent(input$btn_load_example_unpaired_multi, {
    load_example_benchmark("multi_unpaired")
  })
  observeEvent(input$btn_run_uni_eval_atac, {
    run_single_cell_evaluation()
  })
  observeEvent(input$btn_inspect_uni_data_side, { nav_select("datahub_subtabs", "Dataset properties") })
  observeEvent(input$btn_inspect_uni_data_side_atac, { nav_select("datahub_subtabs", "Dataset properties") })
  observeEvent(input$btn_inspect_multi_data_side, { nav_select("datahub_subtabs", "Dataset properties") })
  observeEvent(input$btn_load_demo, { execute_load_demo_benchmark() })
  observeEvent(input$btn_empty_load_demo, { execute_load_demo_benchmark() })
  observeEvent(input$btn_empty_load_demo_view, { execute_load_demo_benchmark() })
  observeEvent(input$btn_empty_load_demo_props, { execute_load_demo_benchmark() })
  observeEvent(input$btn_demo_goto_view, { nav_select("datahub_subtabs", "View Uploaded Data") })
  observeEvent(input$btn_demo_goto_props, { nav_select("datahub_subtabs", "Dataset properties") })
  observeEvent(input$btn_demo_goto_bubble, { nav_select("nav_active", "Comparative Bubble Matrix") })

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # Data Hub: Standalone HPC / CLI Script Generator
  # ----------------------------------------------------------------------------
  generate_hpc_script <- function(mode = "unimodal") {
    date_str <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    is_uni <- (mode == "unimodal")
    n_dsets <- if (is_uni) (input$num_uni_datasets %||% 1) else (input$num_multi_datasets %||% 1)
    n_sims  <- if (is_uni) (input$num_uni_sims %||% 1) else (input$num_multi_sims %||% 1)
    is_paired <- identical(input$opt_multi_pairing, "paired")
    n_cores <- if (is_uni) as.integer(input$num_uni_cores %||% default_n_cores) else as.integer(input$num_multi_cores %||% default_n_cores)
    if (is.na(n_cores) || n_cores < 1) n_cores <- 1L
    
    out <- c(
      "#!/usr/bin/env Rscript",
      "# ==============================================================================",
      "# scSimEval: Standalone Benchmarking Script (User-Configured)",
      sprintf("# Generated on: %s", date_str),
      sprintf("# Benchmark Mode: %s", if (is_uni) "Single-Cell (scRNA or scATAC)" else sprintf("Multiomics (scRNA + scATAC) [%s]", if (is_paired) "Paired" else "Unpaired")),
      sprintf("# Total Datasets: %d | Total Simulators: %d", n_dsets, n_sims),
      sprintf("# Parallel Worker Cores / Threads: %d", n_cores),
      "#",
      "# Instructions:",
      "#   1. Ensure required R packages are installed: install.packages(c('Matrix', 'ggplot2'))",
      "#   2. Ensure scSimEval is installed: remotes::install_local('.') or BiocManager::install('scSimEval')",
      "#   3. Run from terminal or cluster job: Rscript <this_script>.R",
      "# ==============================================================================",
      "",
      "# 1. Load Required Libraries",
      "suppressPackageStartupMessages({",
      "  library(scSimEval)",
      "  library(Matrix)",
      "})",
      "",
      "cat('[scSimEval] Initializing user-defined benchmarking pipeline...\n')",
      "",
      "# Helper: Safe Matrix Loader",
      "load_matrix_safe <- function(path, placeholder = 'matrix') {",
      "  if (is.null(path) || is.na(path) || path == '') {",
      "    stop(sprintf('Path for %s is empty. Please specify a valid file path.', placeholder))",
      "  }",
      "  if (!file.exists(path)) {",
      "    warning(sprintf('File \\\"%s\\\" not found. Using small placeholder matrix for pipeline verification.', path))",
      "    set.seed(42)",
      "    m <- Matrix::rsparsematrix(nrow = 50, ncol = 30, density = 0.2)",
      "    rownames(m) <- paste0('Feature_', seq_len(50))",
      "    colnames(m) <- paste0('Cell_', seq_len(30))",
      "    return(abs(m))",
      "  }",
      "  if (endsWith(tolower(path), '.rds')) {",
      "    m <- readRDS(path)",
      "    if (inherits(m, 'SingleCellExperiment')) {",
      "      m <- SummarizedExperiment::assay(m, 'counts')",
      "    } else if (inherits(m, 'Seurat')) {",
      "      m <- tryCatch(m[['RNA']]$counts, error = function(e) tryCatch(Seurat::GetAssayData(m, slot = 'counts'), error = function(e2) m))",
      "    } else if (is.list(m) && !is.data.frame(m) && !inherits(m, 'Matrix')) {",
      "      candidate_names <- c('new_count', 'counts', 'sim_counts', 'simulated', 'sim_rna', 'sim_atac', 'mat')",
      "      found <- intersect(candidate_names, names(m))",
      "      if (length(found) > 0) m <- m[[found[1]]] else {",
      "        for (item in m) if (is.matrix(item) || inherits(item, 'Matrix') || is.data.frame(item)) { m <- item; break }",
      "      }",
      "    }",
      "  } else {",
      "    sep_char <- if (endsWith(tolower(path), '.csv')) ',' else '\t'",
      "    df <- utils::read.table(path, header = TRUE, row.names = 1, sep = sep_char, check.names = FALSE)",
      "    m <- as.matrix(df)",
      "  }",
      "  if (!inherits(m, 'dgCMatrix')) m <- Matrix::Matrix(m, sparse = TRUE)",
      "  return(m)",
      "}",
      "",
      "# Helper: Safe Vector Loader (Cell Types / Batches)",
      "load_vector_safe <- function(path) {",
      "  if (is.null(path) || is.na(path) || path == '' || !file.exists(path)) return(NULL)",
      "  lower_p <- tolower(path)",
      "  if (endsWith(lower_p, '.rds')) {",
      "    v <- readRDS(path)",
      "    if (is.data.frame(v) || is.matrix(v)) v <- v[, 1]",
      "    return(as.character(v))",
      "  }",
      "  if (endsWith(lower_p, '.xlsx') || endsWith(lower_p, '.xls')) {",
      "    if (requireNamespace('readxl', quietly = TRUE)) {",
      "      df <- as.data.frame(readxl::read_excel(path))",
      "      col_pick <- if ('cell_type' %in% tolower(colnames(df))) which(tolower(colnames(df)) == 'cell_type')[1] else",
      "                  if ('batch' %in% tolower(colnames(df))) which(tolower(colnames(df)) == 'batch')[1] else ncol(df)",
      "      return(as.character(df[[col_pick]]))",
      "    }",
      "  }",
      "  if (endsWith(lower_p, '.csv') || endsWith(lower_p, '.tsv') || endsWith(lower_p, '.txt')) {",
      "    sep_char <- if (endsWith(lower_p, '.csv')) ',' else '\t'",
      "    lines <- readLines(path, n = 5)",
      "    has_header <- length(lines) > 0 && grepl('[a-zA-Z]', lines[1]) && !grepl('^[0-9]+$', lines[1])",
      "    df <- tryCatch(utils::read.table(path, header = has_header, sep = sep_char, stringsAsFactors = FALSE), error = function(e) NULL)",
      "    if (!is.null(df) && is.data.frame(df) && ncol(df) >= 1) {",
      "      col_pick <- if ('cell_type' %in% tolower(colnames(df))) which(tolower(colnames(df)) == 'cell_type')[1] else",
      "                  if ('batch' %in% tolower(colnames(df))) which(tolower(colnames(df)) == 'batch')[1] else ncol(df)",
      "      return(as.character(df[[col_pick]]))",
      "    }",
      "    return(as.character(readLines(path)))",
      "  }",
      "  return(as.character(readLines(path)))",
      "}",
      "",
      "# Helper: Standardize Category Names",
      "standardize_benchmark_categories <- function(df) {",
      "  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(df)",
      "  if (exists('standardize_benchmark_categories', where = asNamespace('scSimEval'), mode = 'function')) {",
      "    return(scSimEval::standardize_benchmark_categories(df))",
      "  }",
      "  df",
      "}",
      "",
      "all_benchmark_results <- list()",
      "all_dataset_summaries <- list()",
      ""
    )
    
    if (is_uni) {
      for (d in seq_len(n_dsets)) {
        dname_custom <- input[[paste0("uni_dset_name_", d)]]
        dname <- if (!is.null(dname_custom) && trimws(dname_custom) != "") trimws(dname_custom) else sprintf("Dataset %d", d)
        ref_f <- input[[paste0("file_uni_ref_", d)]]$name
        if (is.null(ref_f) || ref_f == "") ref_f <- sprintf("ref_%s.rds", make.names(dname))
        ct_f  <- input[[paste0("file_uni_celltypes_", d)]]$name
        bt_f  <- input[[paste0("file_uni_batch_", d)]]$name
        
        out <- c(out,
          sprintf("# ------------------------------------------------------------------------------"),
          sprintf("# Dataset %d: %s", d, dname),
          sprintf("# ------------------------------------------------------------------------------"),
          sprintf("cat('[scSimEval] Processing %s [%d/%d]...\n')", dname, d, n_dsets),
          sprintf("dset_name_%d <- '%s'", d, dname),
          sprintf("ref_file_%d  <- '%s'", d, ref_f),
          sprintf("ref_mat_%d   <- load_matrix_safe(ref_file_%d, placeholder = paste0('Reference ', dset_name_%d))", d, d, d)
        )
        
        if (!is.null(ct_f) && ct_f != "") {
          out <- c(out,
            sprintf("ct_file_%d   <- '%s'", d, ct_f),
            sprintf("ct_vec_%d    <- load_vector_safe(ct_file_%d)", d, d)
          )
        } else {
          out <- c(out, sprintf("ct_vec_%d    <- NULL # (No cell types configured)", d))
        }
        
        if (!is.null(bt_f) && bt_f != "") {
          out <- c(out,
            sprintf("bt_file_%d   <- '%s'", d, bt_f),
            sprintf("bt_vec_%d    <- load_vector_safe(bt_file_%d)", d, d)
          )
        } else {
          out <- c(out, sprintf("bt_vec_%d    <- NULL # (No batches configured)", d))
        }
        
        out <- c(out,
          sprintf("all_dataset_summaries[[length(all_dataset_summaries) + 1]] <- extract_dataset_summary("),
          sprintf("  ref_mat_%d, role = 'Biological Reference', method_name = paste0('Empirical Reference (', dset_name_%d, ')'),", d, d),
          sprintf("  modality = 'Single-Cell (Counts)', cell_types = ct_vec_%d, batch_info = bt_vec_%d", d, d),
          ")",
          ""
        )
        
        for (s in seq_len(n_sims)) {
          sname_custom <- input[[paste0("uni_sim_name_", s)]]
          sname <- if (!is.null(sname_custom) && trimws(sname_custom) != "") trimws(sname_custom) else sprintf("Simulator %d", s)
          sim_f <- input[[paste0("file_uni_sim_", s, "_", d)]]$name
          if (is.null(sim_f) || sim_f == "") sim_f <- sprintf("sim_%s_%s.rds", make.names(sname), make.names(dname))
          sim_t <- as.numeric(input[[paste0("uni_sim_time_", s, "_", d)]]) %||% round(20 + s * 10 + d * 5, 1)
          sim_m <- as.numeric(input[[paste0("uni_sim_mem_", s, "_", d)]]) %||% round(500 + s * 150 + d * 50, 0)
          
          out <- c(out,
            sprintf("  # %s on %s", sname, dname),
            sprintf("  sim_name_%d_%d <- '%s'", s, d, sname),
            sprintf("  sim_file_%d_%d <- '%s'", s, d, sim_f),
            sprintf("  sim_time_%d_%d <- %.1f # Runtime (seconds)", s, d, sim_t),
            sprintf("  sim_mem_%d_%d  <- %.0f # Peak RAM (MiB)", s, d, sim_m),
            sprintf("  cat('[scSimEval]   Evaluating ', sim_name_%d_%d, ' on ', dset_name_%d, '...\n')", s, d, d),
            sprintf("  sim_mat_%d_%d  <- load_matrix_safe(sim_file_%d_%d, placeholder = paste0(sim_name_%d_%d, ' matrix'))", s, d, s, d, s, d),
            "",
            sprintf("  all_dataset_summaries[[length(all_dataset_summaries) + 1]] <- extract_dataset_summary("),
            sprintf("    sim_mat_%d_%d, role = 'Simulated', method_name = paste0(sim_name_%d_%d, ' (', dset_name_%d, ')'),", s, d, s, d, d),
            sprintf("    modality = 'Single-Cell (Counts)', cell_types = ct_vec_%d, batch_info = bt_vec_%d", d, d),
            "  )",
            "",
            sprintf("  res_%d_%d <- evaluate_simulation_accuracy(" , s, d),
            sprintf("    ref_data = ref_mat_%d,", d),
            sprintf("    sim_data = sim_mat_%d_%d,", s, d),
            sprintf("    memory_mb = sim_mem_%d_%d,", s, d),
            sprintf("    elapsed_time = sim_time_%d_%d,", s, d),
            "    compute_bivariate = FALSE,",
            sprintf("    threads = %d,", n_cores),
            "    verbose = FALSE",
            "  )",
            sprintf("  tbl_%d_%d <- res_%d_%d$metrics_summary_table", s, d, s, d),
            sprintf("  tbl_%d_%d$Method  <- sim_name_%d_%d", s, d, s, d),
            sprintf("  tbl_%d_%d$Dataset <- dset_name_%d", s, d, d),
            sprintf("  if (!'Score' %%in%% colnames(tbl_%d_%d) && 'Value' %%in%% colnames(tbl_%d_%d)) tbl_%d_%d$Score <- tbl_%d_%d$Value", s, d, s, d, s, d, s, d),
            sprintf("  all_benchmark_results[[paste0(sim_name_%d_%d, '_', dset_name_%d)]] <- standardize_benchmark_categories(tbl_%d_%d)", s, d, d, s, d),
            ""
          )
        }
      }
    } else {
      # Multiomics mode
      for (d in seq_len(n_dsets)) {
        dname_custom <- input[[paste0("multi_dset_name_", d)]]
        dname <- if (!is.null(dname_custom) && trimws(dname_custom) != "") trimws(dname_custom) else sprintf("Dataset %d", d)
        ref_rna_f  <- input[[paste0("file_multi_ref_rna_", d)]]$name
        ref_atac_f <- input[[paste0("file_multi_ref_atac_", d)]]$name
        if (is.null(ref_rna_f) || ref_rna_f == "") ref_rna_f <- sprintf("ref_rna_%s.rds", make.names(dname))
        if (is.null(ref_atac_f) || ref_atac_f == "") ref_atac_f <- sprintf("ref_atac_%s.rds", make.names(dname))
        ct_f <- input[[paste0("file_multi_celltypes_", d)]]$name
        bt_f <- input[[paste0("file_multi_batch_", d)]]$name
        
        out <- c(out,
          sprintf("# ------------------------------------------------------------------------------"),
          sprintf("# Multiomics Dataset %d: %s", d, dname),
          sprintf("# ------------------------------------------------------------------------------"),
          sprintf("cat('[scSimEval] Processing %s [%d/%d]...\n')", dname, d, n_dsets),
          sprintf("dset_name_%d     <- '%s'", d, dname),
          sprintf("ref_rna_file_%d  <- '%s'", d, ref_rna_f),
          sprintf("ref_atac_file_%d <- '%s'", d, ref_atac_f),
          sprintf("ref_rna_%d       <- load_matrix_safe(ref_rna_file_%d, placeholder = paste0('Reference RNA ', dset_name_%d))", d, d, d),
          sprintf("ref_atac_%d      <- load_matrix_safe(ref_atac_file_%d, placeholder = paste0('Reference ATAC ', dset_name_%d))", d, d, d)
        )
        
        if (!is.null(ct_f) && ct_f != "") {
          out <- c(out,
            sprintf("ct_file_%d        <- '%s'", d, ct_f),
            sprintf("ct_vec_%d         <- load_vector_safe(ct_file_%d)", d, d)
          )
        } else {
          out <- c(out, sprintf("ct_vec_%d         <- NULL # (No cell types configured)", d))
        }
        
        if (!is.null(bt_f) && bt_f != "") {
          out <- c(out,
            sprintf("bt_file_%d        <- '%s'", d, bt_f),
            sprintf("bt_vec_%d         <- load_vector_safe(bt_file_%d)", d, d)
          )
        } else {
          out <- c(out, sprintf("bt_vec_%d         <- NULL # (No batches configured)", d))
        }
        
        out <- c(out,
          sprintf("all_dataset_summaries[[length(all_dataset_summaries) + 1]] <- extract_dataset_summary("),
          sprintf("  ref_rna_%d, role = 'Biological Reference (RNA)', method_name = paste0('Empirical Reference (', dset_name_%d, ')'),", d, d),
          sprintf("  modality = 'RNA-seq', cell_types = ct_vec_%d, batch_info = bt_vec_%d", d, d),
          ")",
          sprintf("all_dataset_summaries[[length(all_dataset_summaries) + 1]] <- extract_dataset_summary("),
          sprintf("  ref_atac_%d, role = 'Biological Reference (ATAC)', method_name = paste0('Empirical Reference (', dset_name_%d, ')'),", d, d),
          sprintf("  modality = 'ATAC-seq', cell_types = ct_vec_%d, batch_info = bt_vec_%d", d, d),
          ")",
          ""
        )
        
        for (s in seq_len(n_sims)) {
          sname_custom <- input[[paste0("multi_sim_name_", s)]]
          sname <- if (!is.null(sname_custom) && trimws(sname_custom) != "") trimws(sname_custom) else sprintf("Simulator %d", s)
          sim_rna_f  <- input[[paste0("file_multi_sim_rna_", s, "_", d)]]$name
          sim_atac_f <- input[[paste0("file_multi_sim_atac_", s, "_", d)]]$name
          if (is.null(sim_rna_f) || sim_rna_f == "") sim_rna_f <- sprintf("sim_rna_%s_%s.rds", make.names(sname), make.names(dname))
          if (is.null(sim_atac_f) || sim_atac_f == "") sim_atac_f <- sprintf("sim_atac_%s_%s.rds", make.names(sname), make.names(dname))
          sim_t <- as.numeric(input[[paste0("multi_sim_time_", s, "_", d)]]) %||% round(45 + s * 15 + d * 10, 1)
          sim_m <- as.numeric(input[[paste0("multi_sim_mem_", s, "_", d)]]) %||% round(1200 + s * 250 + d * 100, 0)
          
          out <- c(out,
            sprintf("  # %s on %s", sname, dname),
            sprintf("  sim_name_%d_%d      <- '%s'", s, d, sname),
            sprintf("  sim_rna_file_%d_%d  <- '%s'", s, d, sim_rna_f),
            sprintf("  sim_atac_file_%d_%d <- '%s'", s, d, sim_atac_f),
            sprintf("  sim_time_%d_%d      <- %.1f # Runtime (seconds)", s, d, sim_t),
            sprintf("  sim_mem_%d_%d       <- %.0f # Peak RAM (MiB)", s, d, sim_m),
            sprintf("  cat('[scSimEval]   Evaluating ', sim_name_%d_%d, ' on ', dset_name_%d, '...\n')", s, d, d),
            sprintf("  sim_rna_%d_%d   <- load_matrix_safe(sim_rna_file_%d_%d, placeholder = paste0(sim_name_%d_%d, ' RNA'))", s, d, s, d, s, d),
            sprintf("  sim_atac_%d_%d  <- load_matrix_safe(sim_atac_file_%d_%d, placeholder = paste0(sim_name_%d_%d, ' ATAC'))", s, d, s, d, s, d),
            "",
            sprintf("  all_dataset_summaries[[length(all_dataset_summaries) + 1]] <- extract_dataset_summary("),
            sprintf("    sim_rna_%d_%d, role = 'Simulated (RNA)', method_name = paste0(sim_name_%d_%d, ' (', dset_name_%d, ')'),", s, d, s, d, d),
            sprintf("    modality = 'RNA-seq', cell_types = ct_vec_%d, batch_info = bt_vec_%d", d, d),
            "  )",
            sprintf("  all_dataset_summaries[[length(all_dataset_summaries) + 1]] <- extract_dataset_summary("),
            sprintf("    sim_atac_%d_%d, role = 'Simulated (ATAC)', method_name = paste0(sim_name_%d_%d, ' (', dset_name_%d, ')'),", s, d, s, d, d),
            sprintf("    modality = 'ATAC-seq', cell_types = ct_vec_%d, batch_info = bt_vec_%d", d, d),
            "  )",
            "",
            sprintf("  res_%d_%d <- evaluate_multiomics_accuracy(", s, d),
            sprintf("    ref_multi    = list(RNA = ref_rna_%d, ATAC = ref_atac_%d),", d, d),
            sprintf("    sim_multi    = list(RNA = sim_rna_%d_%d, ATAC = sim_atac_%d_%d),", s, d, s, d),
            sprintf("    is_paired    = %s,", if (is_paired) "TRUE" else "FALSE"),
            sprintf("    cell_types   = ct_vec_%d,", d),
            sprintf("    batch_info   = bt_vec_%d,", d),
            sprintf("    memory_mb    = sim_mem_%d_%d,", s, d),
            sprintf("    elapsed_time = sim_time_%d_%d,", s, d),
            sprintf("    threads      = %d,", n_cores),
            "    verbose      = FALSE",
            "  )",
            sprintf("  tbl_%d_%d <- if (!is.null(res_%d_%d$metrics_summary_table)) res_%d_%d$metrics_summary_table else res_%d_%d$benchmark_summary_table", s, d, s, d, s, d),
            sprintf("  tbl_%d_%d$Method  <- sim_name_%d_%d", s, d, s, d),
            sprintf("  tbl_%d_%d$Dataset <- dset_name_%d", s, d, d),
            sprintf("  if (!'Score' %%in%% colnames(tbl_%d_%d) && 'Value' %%in%% colnames(tbl_%d_%d)) tbl_%d_%d$Score <- tbl_%d_%d$Value", s, d, s, d, s, d, s, d),
            sprintf("  all_benchmark_results[[paste0(sim_name_%d_%d, '_', dset_name_%d)]] <- standardize_benchmark_categories(tbl_%d_%d)", s, d, d, s, d),
            ""
          )
        }
      }
    }
    
    out <- c(out,
      "# ==============================================================================",
      "# Consolidate & Standardize Benchmark Results",
      "# ==============================================================================",
      "cat('[scSimEval] Consolidating results across methods...\n')",
      "final_benchmark_df <- do.call(rbind, all_benchmark_results)",
      "final_benchmark_df <- standardize_benchmark_categories(final_benchmark_df)",
      "final_summary_df   <- do.call(rbind, all_dataset_summaries)",
      "",
      "# ==============================================================================",
      "# Save Benchmark Data Artifacts",
      "# ==============================================================================",
      "output_rds <- 'scSimEval_benchmark_results.rds'",
      "output_csv <- 'scSimEval_benchmark_results.csv'",
      "output_sum <- 'scSimEval_dataset_summary.csv'",
      "output_pdf <- 'scSimEval_benchmark_bubble_matrix.pdf'",
      "",
      "saveRDS(final_benchmark_df, output_rds)",
      "utils::write.csv(final_benchmark_df, output_csv, row.names = FALSE)",
      "utils::write.csv(final_summary_df, output_sum, row.names = FALSE)",
      "cat('[scSimEval] Results saved to:', output_rds, 'and', output_csv, '\n')",
      "",
      "# ==============================================================================",
      "# Generate Publication Bubble Matrix Figure",
      "# ==============================================================================",
      "if (requireNamespace('ggplot2', quietly = TRUE)) {",
      "  cat('[scSimEval] Generating publication-grade bubble matrix figure...\n')",
      "  p_bubble <- plot_benchmark_bubble_matrix(final_benchmark_df, base_size = 12)",
      "  ggplot2::ggsave(output_pdf, plot = p_bubble, width = 18, height = 9, dpi = 300)",
      "  cat('[scSimEval] Figure exported to:', output_pdf, '\n')",
      "}",
      "",
      "cat('[scSimEval] Benchmarking pipeline completed successfully!\n')",
      "",
      "# ==============================================================================",
      "# ==============================================================================",
      "# FILE PATH CONFIGURATION & HPC CLUSTER RUN GUIDE",
      "# ==============================================================================",
      "# Instructions for Non-Coding Users & HPC / Cluster Deployment:",
      "#   1. Copy this R script and all your data files into your working directory on",
      "#      your HPC cluster, cloud instance, or computing server.",
      "#   2. Ensure all input file variables (ref_*_file, sim_*_file, ct_file, bt_file)",
      "#      above point to your exact file locations on your HPC system.",
      "#   3. Run the script directly from your terminal:",
      "#        Rscript <this_script_name>.R",
      "#      Or submit via SLURM:",
      "#        sbatch --cpus-per-task=8 --mem=32G --time=04:00:00 --wrap=\"Rscript <this_script_name>.R\"",
      "# ==============================================================================",
      ""
    )
    
    return(paste(out, collapse = "\n"))
  }

  output$dl_hpc_r_script_uni <- downloadHandler(
    filename = function() { paste0("run_scSimEval_single_cell_", format(Sys.Date(), "%Y%m%d"), ".R") },
    content  = function(file) { writeLines(generate_hpc_script("unimodal"), file) }
  )
  output$dl_hpc_r_script_uni_side <- downloadHandler(
    filename = function() { paste0("run_scSimEval_single_cell_", format(Sys.Date(), "%Y%m%d"), ".R") },
    content  = function(file) { writeLines(generate_hpc_script("unimodal"), file) }
  )
  output$dl_hpc_r_script_multi <- downloadHandler(
    filename = function() { paste0("run_scSimEval_multiomics_", format(Sys.Date(), "%Y%m%d"), ".R") },
    content  = function(file) { writeLines(generate_hpc_script("multiomics"), file) }
  )
  output$dl_hpc_r_script_multi_side <- downloadHandler(
    filename = function() { paste0("run_scSimEval_multiomics_", format(Sys.Date(), "%Y%m%d"), ".R") },
    content  = function(file) { writeLines(generate_hpc_script("multiomics"), file) }
  )

  # ----------------------------------------------------------------------------
  # Data Hub: Mode 2 Single-Cell Server Logic (Common Reference Datasets 1 to N)
  # ----------------------------------------------------------------------------
  
  # Reactive Parallel Cores Badges
  output$ui_uni_cores_badge <- renderUI({
    n_cores <- as.integer(input$num_uni_cores %||% default_n_cores)
    if (is.na(n_cores) || n_cores < 1) n_cores <- 1L
    span(
      class = "badge px-3 py-2 border shadow-sm",
      style = "background: #F0FDF4; color: #166534; border-color: #BBF7D0 !important; font-size: 0.84rem;",
      icon("microchip", class = "me-1 text-success"),
      sprintf("Parallel: %d CPU Core%s", n_cores, if (n_cores > 1) "s" else "")
    )
  })
  output$ui_uni_cores_badge_bottom <- renderUI({
    n_cores <- as.integer(input$num_uni_cores %||% default_n_cores)
    if (is.na(n_cores) || n_cores < 1) n_cores <- 1L
    span(
      class = "badge px-3 py-2 border shadow-sm",
      style = "background: #F0FDF4; color: #166534; border-color: #BBF7D0 !important; font-size: 0.84rem;",
      icon("microchip", class = "me-1 text-success"),
      sprintf("Parallel: %d CPU Core%s", n_cores, if (n_cores > 1) "s" else "")
    )
  })

  output$ui_multi_cores_badge <- renderUI({
    n_cores <- as.integer(input$num_multi_cores %||% default_n_cores)
    if (is.na(n_cores) || n_cores < 1) n_cores <- 1L
    span(
      class = "badge px-3 py-2 border shadow-sm",
      style = "background: #EEF2FF; color: #3730A3; border-color: #C7D2FE !important; font-size: 0.84rem;",
      icon("microchip", class = "me-1 text-indigo"),
      sprintf("Parallel: %d CPU Core%s", n_cores, if (n_cores > 1) "s" else "")
    )
  })
  output$ui_multi_cores_badge_bottom <- renderUI({
    n_cores <- as.integer(input$num_multi_cores %||% default_n_cores)
    if (is.na(n_cores) || n_cores < 1) n_cores <- 1L
    span(
      class = "badge px-3 py-2 border shadow-sm",
      style = "background: #EEF2FF; color: #3730A3; border-color: #C7D2FE !important; font-size: 0.84rem;",
      icon("microchip", class = "me-1 text-indigo"),
      sprintf("Parallel: %d CPU Core%s", n_cores, if (n_cores > 1) "s" else "")
    )
  })

  output$ui_uni_ref_cards <- renderUI({
    n_dsets <- input$num_uni_datasets %||% 1
    
    dset_cards <- lapply(seq_len(n_dsets), function(d) {
      div(
        class = "mb-3 p-3 rounded",
        style = "background: #F8FAFC; border: 1px solid #CBD5E1; box-shadow: 0 1px 2px rgba(0,0,0,0.03);",
        div(
          class = "d-flex justify-content-between align-items-center mb-2",
          div(
            span(class = "badge bg-primary me-2", sprintf("Reference %d", d)),
            tags$b(sprintf("Biological Reference Dataset %d", d), style = "color: #1E3A8A; font-size: 0.94rem;"),
            tags$small(class = "text-muted ms-2", "(Common ground truth for all simulators)")
          ),
          uiOutput(paste0("ui_uni_ref_badge_", d))
        ),
        fluidRow(
          column(3, textInput(paste0("uni_dset_name_", d), "Dataset Name:", value = sprintf("Dataset %d", d))),
          column(5, fileInput(paste0("file_uni_ref_", d), "Original Reference Count Matrix (.rds / .csv / .tsv / .txt / .xlsx):", accept = c(".rds", ".csv", ".tsv", ".txt", ".xlsx", ".xls"))),
          column(2, fileInput(paste0("file_uni_celltypes_", d), "Cell Types (Optional):", accept = c(".rds", ".csv", ".tsv", ".txt", ".xlsx", ".xls"))),
          column(2, fileInput(paste0("file_uni_batch_", d), "Batches (Optional):", accept = c(".rds", ".csv", ".tsv", ".txt", ".xlsx", ".xls")))
        )
      )
    })
    tagList(dset_cards)
  })
  
  # Badges for Reference Datasets (1 to 20)
  lapply(1:20, function(d) {
    output[[paste0("ui_uni_ref_badge_", d)]] <- renderUI({
      f <- input[[paste0("file_uni_ref_", d)]]
      if (is.null(f) || nrow(f) == 0) return(NULL)
      span(class = "badge bg-success-subtle text-success border border-success-subtle",
           icon("circle-check", class = "me-1"), tools::file_path_sans_ext(f$name))
    })
  })

  # Dynamic Simulator Cards (1 to M) - Under each simulator, upload simulated files for each Reference Dataset
  output$ui_uni_sim_cards <- renderUI({
    n_sims <- input$num_uni_sims %||% 1
    n_dsets <- input$num_uni_datasets %||% 1
    
    sim_cards <- lapply(seq_len(n_sims), function(s) {
      dataset_rows <- lapply(seq_len(n_dsets), function(d) {
        dname_custom <- input[[paste0("uni_dset_name_", d)]]
        dname_disp <- if (!is.null(dname_custom) && trimws(dname_custom) != "") dname_custom else sprintf("Dataset %d", d)
        sim_def <- get_default_sim_scalability(s, d, "unimodal")
        
        div(
          class = "mb-2 p-3 rounded",
          style = "background: #FFFFFF; border: 1px solid #E2E8F0;",
          div(
            class = "d-flex align-items-center mb-2",
            span(class = "badge bg-light text-primary border me-2", sprintf("Dataset %d", d)),
            tags$b(sprintf("Simulated Output for: %s", dname_disp), style = "color: #0F172A; font-size: 0.90rem;"),
            tags$small(class = "text-muted ms-2", "(specify simulation file, runtime, and Peak RAM (MiB) for this dataset)")
          ),
          fluidRow(
            column(6, fileInput(paste0("file_uni_sim_", s, "_", d), sprintf("Simulated Count Matrix for %s (.rds / .csv / .xlsx):", dname_disp), placeholder = sprintf("Simulated %s matrix...", dname_disp), accept = c(".rds", ".csv", ".tsv", ".txt", ".xlsx", ".xls"))),
            column(3, numericInput(paste0("uni_sim_time_", s, "_", d), "Runtime (s):", value = sim_def$time, min = 0.1, step = 1)),
            column(3, numericInput(paste0("uni_sim_mem_", s, "_", d), "Peak RAM (MiB):", value = sim_def$mem, min = 1, step = 10))
          )
        )
      })
      
      tags$details(
        class = "mb-3",
        style = "background: #F8FAFC; border: 1px solid #94A3B8; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.05);",
        open = (s == 1),
        tags$summary(
          style = "padding: 10px 14px; background: #EEF2F6; border-radius: 8px 8px 0 0; cursor: pointer; border-bottom: 1px solid #CBD5E1;",
          div(
            class = "d-flex justify-content-between align-items-center",
            div(
              tags$b(sprintf("Simulator %d", s), style = "font-size: 0.98rem; color: #1E3A8A;"),
              tags$small(class = "text-muted ms-2", "(click to expand/collapse)")
            ),
            span(class = "badge bg-primary", sprintf("%d Dataset%s", n_dsets, if(n_dsets > 1) "s" else ""))
          )
        ),
        div(
          style = "padding: 14px;",
          div(
            style = "max-width: 320px; margin-bottom: 12px;",
            textInput(paste0("uni_sim_name_", s), "Simulator Name:", value = sprintf("Simulator %d", s))
          ),
          tagList(dataset_rows)
        )
      )
    })
    tagList(sim_cards)
  })

  # Single-Cell Evaluation Core Runner Function
  run_single_cell_evaluation <- function() {
    n_dsets <- input$num_uni_datasets %||% 1
    n_sims <- input$num_uni_sims %||% 1
    
    withProgress(message = "Single-Cell Evaluation", value = 0, {
      tryCatch({
        results_list <- list()
        summary_rows <- list()
        all_sim_mats <- list()
        refs_list <- list()
        first_sim_mat <- NULL
        
        for (d in seq_len(n_dsets)) {
          ref_file <- input[[paste0("file_uni_ref_", d)]]
          has_upload_ref <- !is.null(ref_file) && nrow(ref_file) > 0
          
          dname_custom <- input[[paste0("uni_dset_name_", d)]]
          d_label <- if (!is.null(dname_custom) && trimws(dname_custom) != "") trimws(dname_custom) else sprintf("Dataset %d", d)
          
          if (has_upload_ref) {
            incProgress(0.1 / n_dsets, detail = sprintf("Loading reference for %s...", d_label))
            ref_mat_d <- read_uploaded_matrix(ref_file$datapath, ref_file$name)
            if (inherits(ref_mat_d, "matrix") && !inherits(ref_mat_d, "dgCMatrix")) ref_mat_d <- Matrix::Matrix(ref_mat_d, sparse = TRUE)
          } else if (!is.null(initial_demo$toy_data$ref)) {
            ref_mat_d <- initial_demo$toy_data$ref
            d_label <- if (!is.null(dname_custom) && trimws(dname_custom) != "") trimws(dname_custom) else sprintf("Demo Dataset %d", d)
          } else {
            next
          }
          refs_list[[d_label]] <- ref_mat_d
          
          ct_file <- input[[paste0("file_uni_celltypes_", d)]]
          bt_file <- input[[paste0("file_uni_batch_", d)]]
          ct_vec <- if (!is.null(ct_file) && nrow(ct_file) > 0) read_uploaded_labels(ct_file$datapath, ct_file$name) else (if (!has_upload_ref) initial_demo$toy_data$cell_types else NULL)
          bt_vec <- if (!is.null(bt_file) && nrow(bt_file) > 0) read_uploaded_labels(bt_file$datapath, bt_file$name) else (if (!has_upload_ref) initial_demo$toy_data$batch else NULL)
          
          if (!is.null(ct_vec)) {
            rv$cell_types_map[[d_label]] <- ct_vec
            rv$cell_types_map[[paste0(d_label, "_RNA")]] <- ct_vec
          }
          if (!is.null(bt_vec)) {
            rv$batches_map[[d_label]] <- bt_vec
            rv$batches_map[[paste0(d_label, "_RNA")]] <- bt_vec
          }
          
          summary_rows[[length(summary_rows) + 1]] <- extract_dataset_summary(
            ref_mat_d, role = "Biological Reference", method_name = paste0("Empirical Reference (", d_label, ")"),
            modality = "Single-Cell (Counts)", cell_types = ct_vec, batch_info = bt_vec
          )
          
          for (s in seq_len(n_sims)) {
            sim_name <- input[[paste0("uni_sim_name_", s)]]
            if (is.null(sim_name) || trimws(sim_name) == "") sim_name <- sprintf("Simulator %d", s)
            
            sim_t_val <- input[[paste0("uni_sim_time_", s, "_", d)]]
            sim_time <- if (!is.null(sim_t_val) && !is.na(as.numeric(sim_t_val))) as.numeric(sim_t_val) else 30.0
            
            sim_m_val <- input[[paste0("uni_sim_mem_", s, "_", d)]]
            sim_mem <- if (!is.null(sim_m_val) && !is.na(as.numeric(sim_m_val))) as.numeric(sim_m_val) else 800.0
            
            sim_file <- input[[paste0("file_uni_sim_", s, "_", d)]]
            has_upload_sim <- !is.null(sim_file) && nrow(sim_file) > 0
            
            if (has_upload_sim) {
              incProgress(0.7 / (n_sims * n_dsets), detail = sprintf("Evaluating %s on %s...", sim_name, d_label))
              sim_mat_d <- read_uploaded_matrix(sim_file$datapath, sim_file$name)
              if (inherits(sim_mat_d, "matrix") && !inherits(sim_mat_d, "dgCMatrix")) sim_mat_d <- Matrix::Matrix(sim_mat_d, sparse = TRUE)
            } else if (!is.null(initial_demo$toy_data$sim)) {
              incProgress(0.7 / (n_sims * n_dsets), detail = sprintf("Evaluating %s on %s...", sim_name, d_label))
              sim_mat_d <- initial_demo$toy_data$sim
            } else {
              next
            }
            if (is.null(first_sim_mat)) first_sim_mat <- sim_mat_d
            all_sim_mats[[paste0(sim_name, "_", d_label)]] <- sim_mat_d
            
            summary_rows[[length(summary_rows) + 1]] <- extract_dataset_summary(
              sim_mat_d, role = "Simulated", method_name = paste0(sim_name, " (", d_label, ")"),
              modality = "Single-Cell (Counts)", cell_types = ct_vec, batch_info = bt_vec
            )
            
            n_cores_uni <- as.integer(input$num_uni_cores %||% default_n_cores)
            if (is.na(n_cores_uni) || n_cores_uni < 1) n_cores_uni <- 1L
            
            res_sd <- evaluate_simulation_accuracy(
              ref_data = ref_mat_d,
              sim_data = sim_mat_d,
              memory_mb = sim_mem,
              elapsed_time = sim_time,
              compute_bivariate = TRUE,
              threads = n_cores_uni,
              verbose = FALSE
            )
            tbl_sd <- if (!is.null(res_sd$metrics_summary_table)) {
              res_sd$metrics_summary_table
            } else if (!is.null(res_sd$benchmark_summary_table)) {
              res_sd$benchmark_summary_table
            } else {
              res_sd$summary_table
            }
            if (is.null(tbl_sd) || !is.data.frame(tbl_sd)) {
              stop(sprintf("Failed to retrieve single-cell evaluation table for '%s'.", sim_name))
            }
            
            # Align cell types if provided
            ct_sim <- NULL
            ct_ref <- NULL
            if (!is.null(ct_vec) && length(ct_vec) > 0) {
              if (!is.null(names(ct_vec)) && !is.null(colnames(sim_mat_d)) && any(colnames(sim_mat_d) %in% names(ct_vec))) {
                ct_sim <- ct_vec[colnames(sim_mat_d)]
              } else if (length(ct_vec) == ncol(sim_mat_d)) {
                ct_sim <- ct_vec
              }
              if (!is.null(names(ct_vec)) && !is.null(colnames(ref_mat_d)) && any(colnames(ref_mat_d) %in% names(ct_vec))) {
                ct_ref <- ct_vec[colnames(ref_mat_d)]
              } else if (length(ct_vec) == ncol(ref_mat_d)) {
                ct_ref <- ct_vec
              }
            }
            
            # Align batches if provided
            bt_sim <- NULL
            if (!is.null(bt_vec) && length(bt_vec) > 0) {
              if (!is.null(names(bt_vec)) && !is.null(colnames(sim_mat_d)) && any(colnames(sim_mat_d) %in% names(bt_vec))) {
                bt_sim <- bt_vec[colnames(sim_mat_d)]
              } else if (length(bt_vec) == ncol(sim_mat_d)) {
                bt_sim <- bt_vec
              }
            }
            
            # Cellular Structure & Concordance (Category III)
            if (!is.null(ct_sim) && length(ct_sim) > 0 && length(unique(stats::na.omit(ct_sim))) > 1) {
              try({
                clu <- evaluate_clustering_metrics(data = sim_mat_d, cell_types = ct_sim, ref_data = ref_mat_d)
                clu_names <- c("silhouette_sim", "dunn_sim", "davies_bouldin_sim", "calinski_harabasz_sim", "ari", "nmi", "ami", "fmi", "homogeneity", "completeness", "v_measure", "neighborhood_purity")
                clu_rows <- list()
                for (cm in clu_names) {
                  if (!is.null(clu[[cm]])) {
                    val_cm <- as.numeric(clu[[cm]])
                    if (!is.na(val_cm)) {
                      clu_rows[[length(clu_rows) + 1]] <- data.frame(
                        Category = "(III) Cellular Structure & Concordance",
                        Property = "clustering_concordance",
                        Metric = cm,
                        Value = val_cm,
                        stringsAsFactors = FALSE
                      )
                    }
                  }
                }
                if (length(clu_rows) > 0) tbl_sd <- rbind(tbl_sd, do.call(rbind, clu_rows))
              }, silent = TRUE)
            }
            
            # Batch Effects & Confounder Mixing (Category IV)
            if (!is.null(bt_sim) && length(bt_sim) > 0 && length(unique(stats::na.omit(bt_sim))) > 1) {
              try({
                sim_m_num <- as.matrix(sim_mat_d)
                mode(sim_m_num) <- "double"
                bat <- evaluate_batch_metrics(data = sim_m_num, batch_info = bt_sim, cell_types = ct_sim)
                bat_names <- c("batch_silhouette", "shannon_entropy", "pcr_r2", "cms", "isi", "seurat_mixing_metric", "ldf_diff")
                bat_rows <- list()
                for (bm in bat_names) {
                  if (!is.null(bat[[bm]])) {
                    val_bm <- as.numeric(bat[[bm]])
                    if (!is.na(val_bm)) {
                      bat_rows[[length(bat_rows) + 1]] <- data.frame(
                        Category = "(IV) Batch Effects & Confounder Mixing",
                        Property = "batch_mixing",
                        Metric = bm,
                        Value = val_bm,
                        stringsAsFactors = FALSE
                      )
                    }
                  }
                }
                if (length(bat_rows) > 0) tbl_sd <- rbind(tbl_sd, do.call(rbind, bat_rows))
              }, silent = TRUE)
            }
            
            # Biological Signal & Downstream Fidelity (Category V)
            if (!is.null(ct_sim) && !is.null(ct_ref) && length(unique(stats::na.omit(ct_sim))) > 1 && length(unique(stats::na.omit(ct_ref))) > 1) {
              try({
                sig_eval <- evaluate_simbench_signals(ref_mat = ref_mat_d, sim_mat = sim_mat_d, ref_celltypes = ct_ref, sim_celltypes = ct_sim)
                if (is.data.frame(sig_eval) && nrow(sig_eval) > 0) {
                  sig_rows <- list()
                  for (i in seq_len(nrow(sig_eval))) {
                    val_sig <- as.numeric(sig_eval$Absolute_Error[i])
                    if (!is.na(val_sig)) {
                      sig_rows[[length(sig_rows) + 1]] <- data.frame(
                        Category = "(V) Biological Signal & Downstream Fidelity",
                        Property = "simbench_signal",
                        Metric = paste0(sig_eval$Signal_Type[i], "_fidelity"),
                        Value = val_sig,
                        stringsAsFactors = FALSE
                      )
                    }
                  }
                  if (length(sig_rows) > 0) tbl_sd <- rbind(tbl_sd, do.call(rbind, sig_rows))
                }
              }, silent = TRUE)
            }
            
            tbl_sd$Method <- sim_name
            tbl_sd$Dataset <- d_label
            if ("Value" %in% colnames(tbl_sd)) tbl_sd$Score_Raw <- tbl_sd$Value
            if (!"Score" %in% colnames(tbl_sd) && "Value" %in% colnames(tbl_sd)) tbl_sd$Score <- tbl_sd$Value
            if (!"Score_Norm" %in% colnames(tbl_sd) && "Score" %in% colnames(tbl_sd)) tbl_sd$Score_Norm <- tbl_sd$Score
            results_list[[paste0(sim_name, "_", d_label)]] <- standardize_benchmark_categories(tbl_sd)
            rm(sim_mat_d); gc(verbose = FALSE)
          }
        }
        
        if (length(results_list) == 0) stop("Please upload reference matrix and at least one simulated matrix.")
        
        combined_df <- do.call(rbind, results_list)
        combined_df <- standardize_benchmark_categories(combined_df)
        
        if (isTRUE(input$chk_append_uni) && !is.null(rv$benchmark_df)) {
          existing_clean <- rv$benchmark_df[!rv$benchmark_df$Method %in% unique(combined_df$Method), , drop = FALSE]
          rv$benchmark_df <- rbind(existing_clean, combined_df[, intersect(colnames(existing_clean), colnames(combined_df))])
        } else {
          rv$benchmark_df <- combined_df
        }
        
        new_summary_df <- do.call(rbind, summary_rows)
        if (isTRUE(input$chk_append_uni) && !is.null(rv$dataset_summary_df)) {
          rv$dataset_summary_df <- rbind(rv$dataset_summary_df, new_summary_df)
        } else {
          rv$dataset_summary_df <- new_summary_df
        }
        
        if (!is.null(refs_list) && length(refs_list) > 0) rv$toy_ref <- refs_list[[1]]
        if (length(all_sim_mats) > 0) rv$sim_matrices <- all_sim_mats
        rv$ref_matrices <- refs_list
        if (!is.null(ct_vec)) rv$cell_types <- ct_vec
        if (!is.null(bt_vec)) rv$batch <- bt_vec
        
        # Consolidate reference and simulated matrices for inspection
        all_combined_mats <- list()
        for (r_name in names(refs_list)) {
          all_combined_mats[[paste0("[Reference] ", r_name)]] <- refs_list[[r_name]]
        }
        for (s_name in names(all_sim_mats)) {
          all_combined_mats[[paste0("[Simulated] ", s_name)]] <- all_sim_mats[[s_name]]
        }
        rv$all_matrices <- all_combined_mats
        
        mode_data$unimodal <- list(
          benchmark_df       = rv$benchmark_df,
          methods            = unique(rv$benchmark_df$Method),
          all_matrices       = all_combined_mats,
          ref_matrices       = refs_list,
          sim_matrices       = all_sim_mats,
          dataset_summary_df = rv$dataset_summary_df,
          toy_ref            = rv$toy_ref,
          toy_sim            = rv$toy_sim,
          cell_types         = rv$cell_types,
          batch              = rv$batch,
          cell_types_map     = rv$cell_types_map,
          batches_map        = rv$batches_map,
          source_name        = "User Evaluated Single-Cell Benchmarks"
        )
        sync_active_mode_data("unimodal")
        
        sc_log("SUCCESS", paste0("Single-Cell Evaluation completed successfully: ", nrow(rv$benchmark_df), " metric rows generated across ", length(unique(rv$benchmark_df$Method)), " method(s)."), "Evaluation")
        showNotification("Single-Cell Evaluation completed successfully!", type = "message", duration = 5)
      }, error = function(e) {
        sc_log("ERROR", paste("Single-Cell Evaluation error:", e$message), "Evaluation")
        showNotification(paste("Error during evaluation:", e$message), type = "error", duration = 8)
      })
    })
  }

  observeEvent(input$btn_run_uni_eval, run_single_cell_evaluation())
  observeEvent(input$btn_run_uni_eval_main, run_single_cell_evaluation())

  # ----------------------------------------------------------------------------
  # Data Hub: Mode 3 Multiomics Server Logic (Common Reference Datasets 1 to N)
  # ----------------------------------------------------------------------------
  output$ui_multi_ref_cards <- renderUI({
    n_dsets <- input$num_multi_datasets %||% 1
    is_unpaired <- identical(input$opt_multi_pairing, "unpaired")
    
    dset_cards <- lapply(seq_len(n_dsets), function(d) {
      div(
        class = "mb-3 p-3 rounded",
        style = "background: #F8FAFC; border: 1px solid #CBD5E1; box-shadow: 0 1px 2px rgba(0,0,0,0.03);",
        div(
          class = "d-flex justify-content-between align-items-center mb-2",
          div(
            span(class = "badge text-white me-2", style = "background: #4F46E5;", sprintf("Reference %d", d)),
            tags$b(sprintf("Biological Reference Multiomics Dataset %d", d), style = "color: #312E81; font-size: 0.94rem;"),
            tags$small(class = "text-muted ms-2", if (is_unpaired) "(Unpaired: independent cells for RNA and ATAC)" else "(Paired single-cell co-assay across all multiomics simulators)")
          ),
          div(uiOutput(paste0("ui_multi_ref_rna_badge_", d)), uiOutput(paste0("ui_multi_ref_atac_badge_", d)))
        ),
        if (is_unpaired) {
          tagList(
            fluidRow(
              column(2, textInput(paste0("multi_dset_name_", d), "Dataset Name:", value = sprintf("Dataset %d", d))),
              column(5, fileInput(paste0("file_multi_ref_rna_", d), "Reference RNA Matrix (.rds / .csv):", accept = c(".rds", ".csv", ".tsv", ".txt"))),
              column(5, fileInput(paste0("file_multi_ref_atac_", d), "Reference ATAC Matrix (.rds / .csv):", accept = c(".rds", ".csv", ".tsv", ".txt")))
            ),
            div(
              style = "background: #EEF2FF; border: 1px dashed #A5B4FC; border-radius: 6px; padding: 10px; margin-top: 6px;",
              div(style = "font-size: 0.78rem; font-weight: 700; color: #3730A3; margin-bottom: 6px;",
                  icon("layer-group", class = "me-1"), "Unpaired Annotations (Independent RNA & ATAC Cells):"),
              fluidRow(
                column(3, fileInput(paste0("file_multi_celltypes_rna_", d), "RNA Cell Types (.rds / .csv):", accept = c(".rds", ".csv", ".tsv", ".txt"))),
                column(3, fileInput(paste0("file_multi_celltypes_atac_", d), "ATAC Cell Types (.rds / .csv):", accept = c(".rds", ".csv", ".tsv", ".txt"))),
                column(3, fileInput(paste0("file_multi_batch_rna_", d), "RNA Batches (.rds / .csv):", accept = c(".rds", ".csv", ".tsv", ".txt"))),
                column(3, fileInput(paste0("file_multi_batch_atac_", d), "ATAC Batches (.rds / .csv):", accept = c(".rds", ".csv", ".tsv", ".txt")))
              )
            )
          )
        } else {
          fluidRow(
            column(2, textInput(paste0("multi_dset_name_", d), "Dataset Name:", value = sprintf("Dataset %d", d))),
            column(3, fileInput(paste0("file_multi_ref_rna_", d), "Reference RNA Matrix (.rds / .csv):", accept = c(".rds", ".csv", ".tsv", ".txt"))),
            column(3, fileInput(paste0("file_multi_ref_atac_", d), "Reference ATAC Matrix (.rds / .csv):", accept = c(".rds", ".csv", ".tsv", ".txt"))),
            column(2, fileInput(paste0("file_multi_celltypes_", d), "Cell Types (Opt):", accept = c(".rds", ".csv", ".tsv", ".txt"))),
            column(2, fileInput(paste0("file_multi_batch_", d), "Batches (Opt):", accept = c(".rds", ".csv", ".tsv", ".txt")))
          )
        }
      )
    })
    tagList(dset_cards)
  })
  
  # Badges for Multiomics Reference Datasets
  lapply(1:20, function(d) {
    output[[paste0("ui_multi_ref_rna_badge_", d)]] <- renderUI({
      f <- input[[paste0("file_multi_ref_rna_", d)]]
      if (is.null(f) || nrow(f) == 0) return(NULL)
      span(class = "badge bg-primary text-white me-1", icon("dna", class = "me-1"), tools::file_path_sans_ext(f$name))
    })
    output[[paste0("ui_multi_ref_atac_badge_", d)]] <- renderUI({
      f <- input[[paste0("file_multi_ref_atac_", d)]]
      if (is.null(f) || nrow(f) == 0) return(NULL)
      span(class = "badge bg-info text-white", icon("chart-simple", class = "me-1"), tools::file_path_sans_ext(f$name))
    })
  })

  # Dynamic Multiomics Simulator Cards (1 to M)
  output$ui_multi_sim_cards <- renderUI({
    n_sims <- input$num_multi_sims %||% 1
    n_dsets <- input$num_multi_datasets %||% 1
    is_unpaired <- identical(input$opt_multi_pairing, "unpaired")
    
    sim_cards <- lapply(seq_len(n_sims), function(s) {
      dataset_rows <- lapply(seq_len(n_dsets), function(d) {
        dname_custom <- input[[paste0("multi_dset_name_", d)]]
        dname_disp <- if (!is.null(dname_custom) && trimws(dname_custom) != "") dname_custom else sprintf("Dataset %d", d)
        sim_def <- get_default_sim_scalability(s, d, "multiomics")
        
        div(
          class = "mb-2 p-3 rounded",
          style = "background: #FFFFFF; border: 1px solid #E2E8F0;",
          div(
            class = "d-flex align-items-center mb-2",
            span(class = "badge text-white me-2", style = "background: #4F46E5;", sprintf("Dataset %d", d)),
            tags$b(sprintf("Simulated Multiomics Data for: %s", dname_disp), style = "color: #312E81; font-size: 0.90rem;"),
            tags$small(class = "text-muted ms-2", "(specify simulated RNA & ATAC, runtime, and Peak RAM (MiB) for this dataset)")
          ),
          if (is_unpaired) {
            tagList(
              fluidRow(
                column(6, fileInput(paste0("file_multi_sim_rna_", s, "_", d), sprintf("Simulated RNA for %s (.rds / .csv):", dname_disp), accept = c(".rds", ".csv", ".tsv", ".txt"))),
                column(6, fileInput(paste0("file_multi_sim_atac_", s, "_", d), sprintf("Simulated ATAC for %s (.rds / .csv):", dname_disp), accept = c(".rds", ".csv", ".tsv", ".txt")))
              ),
              div(
                style = "background: #F8FAFC; border: 1px dashed #CBD5E1; border-radius: 6px; padding: 8px 12px; margin-top: 4px;",
                div(style = "font-size: 0.76rem; font-weight: 700; color: #475569; margin-bottom: 4px;",
                    icon("stopwatch", class = "me-1 text-primary"), "Dual / Flexible Scalability Metrics (Unpaired Simulator):"),
                fluidRow(
                  column(3, numericInput(paste0("multi_sim_time_rna_", s, "_", d), "RNA Runtime (s):", value = round(sim_def$time * 0.45, 1), min = 0.1, step = 1)),
                  column(3, numericInput(paste0("multi_sim_mem_rna_", s, "_", d), "RNA Peak RAM (MiB):", value = round(sim_def$mem * 0.85, 0), min = 1, step = 10)),
                  column(3, numericInput(paste0("multi_sim_time_atac_", s, "_", d), "ATAC Runtime (s):", value = round(sim_def$time * 0.55, 1), min = 0.1, step = 1)),
                  column(3, numericInput(paste0("multi_sim_mem_atac_", s, "_", d), "ATAC Peak RAM (MiB):", value = round(sim_def$mem * 0.95, 0), min = 1, step = 10))
                )
              )
            )
          } else {
            fluidRow(
              column(4, fileInput(paste0("file_multi_sim_rna_", s, "_", d), sprintf("Simulated RNA for %s (.rds / .csv):", dname_disp), accept = c(".rds", ".csv", ".tsv", ".txt"))),
              column(4, fileInput(paste0("file_multi_sim_atac_", s, "_", d), sprintf("Simulated ATAC for %s (.rds / .csv):", dname_disp), accept = c(".rds", ".csv", ".tsv", ".txt"))),
              column(2, numericInput(paste0("multi_sim_time_", s, "_", d), "Runtime (s):", value = sim_def$time, min = 0.1, step = 1)),
              column(2, numericInput(paste0("multi_sim_mem_", s, "_", d), "Peak RAM (MiB):", value = sim_def$mem, min = 1, step = 10))
            )
          }
        )
      })
      
      tags$details(
        class = "mb-3",
        style = "background: #F8FAFC; border: 1px solid #94A3B8; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.05);",
        open = (s == 1),
        tags$summary(
          style = "padding: 10px 14px; background: #EEF2FF; border-radius: 8px 8px 0 0; cursor: pointer; border-bottom: 1px solid #C7D2FE;",
          div(
            class = "d-flex justify-content-between align-items-center",
            div(
              tags$b(sprintf("Simulator %d", s), style = "font-size: 0.98rem; color: #312E81;"),
              tags$small(class = "text-muted ms-2", "(click to expand/collapse)")
            ),
            span(class = "badge text-white", style = "background: #4F46E5;", sprintf("%d Dataset%s", n_dsets, if(n_dsets > 1) "s" else ""))
          )
        ),
        div(
          style = "padding: 14px;",
          div(
            style = "max-width: 320px; margin-bottom: 12px;",
            textInput(paste0("multi_sim_name_", s), "Simulator Name:", value = sprintf("Simulator %d", s))
          ),
          tagList(dataset_rows)
        )
      )
    })
    tagList(sim_cards)
  })

  # Multiomics Evaluation Core Runner Function
  run_multiomics_evaluation <- function() {
    n_dsets <- input$num_multi_datasets %||% 1
    n_sims <- input$num_multi_sims %||% 1
    is_paired <- identical(input$opt_multi_pairing, "paired")
    
    withProgress(message = "Multiomics Evaluation", value = 0, {
      tryCatch({
        results_list <- list()
        summary_rows <- list()
        all_sim_mats <- list()
        refs_list <- list()
        first_sim_mat <- NULL
        
        for (d in seq_len(n_dsets)) {
          ref_rna_file  <- input[[paste0("file_multi_ref_rna_", d)]]
          ref_atac_file <- input[[paste0("file_multi_ref_atac_", d)]]
          if (is.null(ref_rna_file) || is.null(ref_atac_file)) next
          
          dname_custom <- input[[paste0("multi_dset_name_", d)]]
          d_label <- if (!is.null(dname_custom) && trimws(dname_custom) != "") dname_custom else sprintf("Dataset %d", d)
          
          incProgress(0.1 / n_dsets, detail = sprintf("Loading reference RNA & ATAC for %s...", d_label))
          ref_rna_d <- read_uploaded_matrix(ref_rna_file$datapath, ref_rna_file$name)
          ref_atac_d <- read_uploaded_matrix(ref_atac_file$datapath, ref_atac_file$name)
          if (inherits(ref_rna_d, "matrix") && !inherits(ref_rna_d, "dgCMatrix")) ref_rna_d <- Matrix::Matrix(ref_rna_d, sparse = TRUE)
          if (inherits(ref_atac_d, "matrix") && !inherits(ref_atac_d, "dgCMatrix")) ref_atac_d <- Matrix::Matrix(ref_atac_d, sparse = TRUE)
          refs_list[[paste0(d_label, "_RNA")]] <- ref_rna_d
          refs_list[[paste0(d_label, "_ATAC")]] <- ref_atac_d
          
          if (is_paired) {
            ct_file <- input[[paste0("file_multi_celltypes_", d)]]
            bt_file <- input[[paste0("file_multi_batch_", d)]]
            ct_vec <- if (!is.null(ct_file) && nrow(ct_file) > 0) read_uploaded_labels(ct_file$datapath, ct_file$name) else NULL
            bt_vec <- if (!is.null(bt_file) && nrow(bt_file) > 0) read_uploaded_labels(bt_file$datapath, bt_file$name) else NULL
            eval_ct <- ct_vec
            eval_bt <- bt_vec
            ct_r_disp <- ct_vec; ct_a_disp <- ct_vec
            bt_r_disp <- bt_vec; bt_a_disp <- bt_vec
          } else {
            ct_r_f <- input[[paste0("file_multi_celltypes_rna_", d)]] %||% input[[paste0("file_multi_celltypes_", d)]]
            ct_a_f <- input[[paste0("file_multi_celltypes_atac_", d)]]
            bt_r_f <- input[[paste0("file_multi_batch_rna_", d)]] %||% input[[paste0("file_multi_batch_", d)]]
            bt_a_f <- input[[paste0("file_multi_batch_atac_", d)]]
            
            ct_r_disp <- if (!is.null(ct_r_f) && nrow(ct_r_f) > 0) read_uploaded_labels(ct_r_f$datapath, ct_r_f$name) else NULL
            ct_a_disp <- if (!is.null(ct_a_f) && nrow(ct_a_f) > 0) read_uploaded_labels(ct_a_f$datapath, ct_a_f$name) else NULL
            bt_r_disp <- if (!is.null(bt_r_f) && nrow(bt_r_f) > 0) read_uploaded_labels(bt_r_f$datapath, bt_r_f$name) else NULL
            bt_a_disp <- if (!is.null(bt_a_f) && nrow(bt_a_f) > 0) read_uploaded_labels(bt_a_f$datapath, bt_a_f$name) else NULL
            
            eval_ct <- if (!is.null(ct_r_disp) || !is.null(ct_a_disp)) list(RNA = ct_r_disp, ATAC = ct_a_disp) else NULL
            eval_bt <- if (!is.null(bt_r_disp) || !is.null(bt_a_disp)) list(RNA = bt_r_disp, ATAC = bt_a_disp) else NULL
          }
          
          if (!is.null(ct_r_disp)) {
            rv$cell_types_map[[paste0(d_label, "_RNA")]] <- ct_r_disp
            rv$cell_types_map[[d_label]] <- ct_r_disp
          }
          if (!is.null(ct_a_disp)) {
            rv$cell_types_map[[paste0(d_label, "_ATAC")]] <- ct_a_disp
          }
          if (!is.null(bt_r_disp)) {
            rv$batches_map[[paste0(d_label, "_RNA")]] <- bt_r_disp
            rv$batches_map[[d_label]] <- bt_r_disp
          }
          if (!is.null(bt_a_disp)) {
            rv$batches_map[[paste0(d_label, "_ATAC")]] <- bt_a_disp
          }
          if (d == 1) {
            rv$cell_types <- if (!is.null(ct_r_disp)) ct_r_disp else ct_a_disp
            rv$batch <- if (!is.null(bt_r_disp)) bt_r_disp else bt_a_disp
          }
          
          summary_rows[[length(summary_rows) + 1]] <- extract_dataset_summary(
            ref_rna_d, role = "Biological Reference (RNA)", method_name = paste0("Empirical Reference (", d_label, ")"),
            modality = "RNA-seq", cell_types = ct_r_disp, batch_info = bt_r_disp
          )
          summary_rows[[length(summary_rows) + 1]] <- extract_dataset_summary(
            ref_atac_d, role = "Biological Reference (ATAC)", method_name = paste0("Empirical Reference (", d_label, ")"),
            modality = "ATAC-seq", cell_types = ct_a_disp, batch_info = bt_a_disp
          )
          
          for (s in seq_len(n_sims)) {
            sim_name <- input[[paste0("multi_sim_name_", s)]]
            if (is.null(sim_name) || trimws(sim_name) == "") sim_name <- sprintf("Simulator %d", s)
            
            sim_rna_file  <- input[[paste0("file_multi_sim_rna_", s, "_", d)]]
            sim_atac_file <- input[[paste0("file_multi_sim_atac_", s, "_", d)]]
            if (is.null(sim_rna_file) || is.null(sim_atac_file)) next
            
            t_r <- input[[paste0("multi_sim_time_rna_", s, "_", d)]]
            t_a <- input[[paste0("multi_sim_time_atac_", s, "_", d)]]
            m_r <- input[[paste0("multi_sim_mem_rna_", s, "_", d)]]
            m_a <- input[[paste0("multi_sim_mem_atac_", s, "_", d)]]
            
            sim_time <- if (!is_paired && !is.null(t_r) && !is.null(t_a)) {
              as.numeric(t_r) + as.numeric(t_a)
            } else {
              as.numeric(input[[paste0("multi_sim_time_", s, "_", d)]]) %||% 60.0
            }
            sim_mem <- if (!is_paired && !is.null(m_r) && !is.null(m_a)) {
              max(as.numeric(m_r), as.numeric(m_a))
            } else {
              as.numeric(input[[paste0("multi_sim_mem_", s, "_", d)]]) %||% 1500.0
            }
            
            incProgress(0.7 / (n_sims * n_dsets), detail = sprintf("Evaluating %s on %s...", sim_name, d_label))
            sim_rna_d  <- read_uploaded_matrix(sim_rna_file$datapath, sim_rna_file$name)
            sim_atac_d <- read_uploaded_matrix(sim_atac_file$datapath, sim_atac_file$name)
            if (inherits(sim_rna_d, "matrix") && !inherits(sim_rna_d, "dgCMatrix")) sim_rna_d <- Matrix::Matrix(sim_rna_d, sparse = TRUE)
            if (inherits(sim_atac_d, "matrix") && !inherits(sim_atac_d, "dgCMatrix")) sim_atac_d <- Matrix::Matrix(sim_atac_d, sparse = TRUE)
            if (is.null(first_sim_mat)) first_sim_mat <- sim_rna_d
            all_sim_mats[[paste0(sim_name, "_RNA_", d_label)]] <- sim_rna_d
            all_sim_mats[[paste0(sim_name, "_ATAC_", d_label)]] <- sim_atac_d
            
            summary_rows[[length(summary_rows) + 1]] <- extract_dataset_summary(
              sim_rna_d, role = "Simulated (RNA)", method_name = paste0(sim_name, " (", d_label, ")"),
              modality = "RNA-seq", cell_types = ct_r_disp, batch_info = bt_r_disp
            )
            summary_rows[[length(summary_rows) + 1]] <- extract_dataset_summary(
              sim_atac_d, role = "Simulated (ATAC)", method_name = paste0(sim_name, " (", d_label, ")"),
              modality = "ATAC-seq", cell_types = ct_a_disp, batch_info = bt_a_disp
            )
            
            n_cores_multi <- as.integer(input$num_multi_cores %||% default_n_cores)
            if (is.na(n_cores_multi) || n_cores_multi < 1) n_cores_multi <- 1L
            
            res_multi_sd <- evaluate_multiomics_accuracy(
              ref_multi = list(RNA = ref_rna_d, ATAC = ref_atac_d),
              sim_multi = list(RNA = sim_rna_d, ATAC = sim_atac_d),
              is_paired = is_paired,
              cell_types = eval_ct,
              batch_info = eval_bt,
              memory_mb = sim_mem,
              elapsed_time = sim_time,
              threads = n_cores_multi,
              verbose = FALSE
            )
            tbl_m_sd <- if (!is.null(res_multi_sd$metrics_summary_table)) {
              res_multi_sd$metrics_summary_table
            } else if (!is.null(res_multi_sd$benchmark_summary_table)) {
              res_multi_sd$benchmark_summary_table
            } else {
              res_multi_sd$summary_table
            }
            if (is.null(tbl_m_sd) || !is.data.frame(tbl_m_sd)) {
              stop(sprintf("Failed to retrieve multiomics evaluation table for '%s'.", sim_name))
            }
            tbl_m_sd$Method <- sim_name
            tbl_m_sd$Dataset <- d_label
            if ("Value" %in% colnames(tbl_m_sd)) tbl_m_sd$Score_Raw <- tbl_m_sd$Value
            if (!"Score" %in% colnames(tbl_m_sd) && "Value" %in% colnames(tbl_m_sd)) tbl_m_sd$Score <- tbl_m_sd$Value
            if (!"Score_Norm" %in% colnames(tbl_m_sd) && "Score" %in% colnames(tbl_m_sd)) tbl_m_sd$Score_Norm <- tbl_m_sd$Score
            results_list[[paste0(sim_name, "_", d_label)]] <- standardize_benchmark_categories(tbl_m_sd)
            rm(sim_rna_d, sim_atac_d); gc(verbose = FALSE)
          }
        }
        
        if (length(results_list) == 0) stop("Please upload reference RNA/ATAC and simulated files.")
        
        combined_df <- do.call(rbind, results_list)
        combined_df <- standardize_benchmark_categories(combined_df)
        
        if (isTRUE(input$chk_append_multi) && !is.null(rv$benchmark_df)) {
          existing_clean <- rv$benchmark_df[!rv$benchmark_df$Method %in% unique(combined_df$Method), , drop = FALSE]
          rv$benchmark_df <- rbind(existing_clean, combined_df[, intersect(colnames(existing_clean), colnames(combined_df))])
        } else {
          rv$benchmark_df <- combined_df
        }
        
        new_summary_df <- do.call(rbind, summary_rows)
        if (isTRUE(input$chk_append_multi) && !is.null(rv$dataset_summary_df)) {
          rv$dataset_summary_df <- rbind(rv$dataset_summary_df, new_summary_df)
        } else {
          rv$dataset_summary_df <- new_summary_df
        }
        
        if (!is.null(refs_list) && length(refs_list) > 0) rv$toy_ref <- refs_list[[1]]
        if (!is.null(first_sim_mat)) rv$toy_sim <- first_sim_mat
        if (length(all_sim_mats) > 0) rv$sim_matrices <- all_sim_mats
        rv$ref_matrices <- refs_list
        
        # Consolidate reference and simulated multiomics matrices for inspection
        all_combined_multi <- list()
        for (r_name in names(refs_list)) {
          all_combined_multi[[paste0("[Reference] ", r_name)]] <- refs_list[[r_name]]
        }
        for (s_name in names(all_sim_mats)) {
          all_combined_multi[[paste0("[Simulated] ", s_name)]] <- all_sim_mats[[s_name]]
        }
        rv$all_matrices <- all_combined_multi
        
        mode_data$multiomics <- list(
          benchmark_df       = rv$benchmark_df,
          methods            = unique(rv$benchmark_df$Method),
          all_matrices       = all_combined_multi,
          ref_matrices       = refs_list,
          sim_matrices       = all_sim_mats,
          dataset_summary_df = rv$dataset_summary_df,
          toy_ref            = rv$toy_ref,
          toy_sim            = rv$toy_sim,
          cell_types         = rv$cell_types,
          batch              = rv$batch,
          cell_types_map     = rv$cell_types_map,
          batches_map        = rv$batches_map,
          source_name        = "User Evaluated Multiomics Benchmarks"
        )
        sync_active_mode_data("multiomics")
        
        sc_log("SUCCESS", paste0("Multiomics Evaluation completed successfully: ", nrow(rv$benchmark_df), " evaluation rows across ", length(unique(rv$benchmark_df$Method)), " method(s)."), "Evaluation")
        showNotification("Multiomics Evaluation completed successfully!", type = "message", duration = 5)
      }, error = function(e) {
        sc_log("ERROR", paste("Multiomics Evaluation error:", e$message), "Evaluation")
        showNotification(paste("Error during multiomics evaluation:", e$message), type = "error", duration = 8)
      })
    })
  }

  observeEvent(input$btn_run_multi_eval, run_multiomics_evaluation())
  observeEvent(input$btn_run_multi_eval_main, run_multiomics_evaluation())
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
        
        if (is.list(obj) && !is.null(obj$dataset_summary)) {
          rv$dataset_summary_df <- obj$dataset_summary
        } else if (is.list(obj) && !is.null(obj$dataset_properties)) {
          rv$dataset_summary_df <- obj$dataset_properties
        } else {
          rv$dataset_summary_df <- NULL
        }
      } else if (ext == "csv") {
        df <- utils::read.csv(input$file_bench_upload$datapath, check.names = FALSE)
        rv$dataset_summary_df <- NULL
      } else {
        stop("Unsupported file type. Please upload .rds or .csv.")
      }
      
      if (!"Method" %in% colnames(df) && "Simulator" %in% colnames(df)) df$Method <- df$Simulator
      if (!"Score" %in% colnames(df) && "Value" %in% colnames(df)) df$Score <- df$Value
      df <- standardize_benchmark_categories(df)
      
      rv$benchmark_df <- df
      rv$methods <- unique(df$Method)
      rv$source_name <- paste0("Uploaded File: ", input$file_bench_upload$name)
      mode_data$upload_bench <- list(
        benchmark_df       = df,
        methods            = rv$methods,
        all_matrices       = list(),
        ref_matrices       = list(),
        sim_matrices       = list(),
        dataset_summary_df = NULL,
        toy_ref            = NULL,
        toy_sim            = NULL,
        cell_types         = NULL,
        batch              = NULL,
        source_name        = rv$source_name
      )
      sync_active_mode_data("upload_bench")
      
      updateCheckboxGroupInput(session, "sel_bubble_methods", choices = rv$methods, selected = rv$methods)
      sc_log("SUCCESS", paste0("Uploaded benchmark results loaded from '", input$file_bench_upload$name, "': ", nrow(df), " records across ", length(rv$methods), " simulators: ", paste(rv$methods, collapse = ", "), "."), "Upload")
      showNotification("Benchmark results loaded successfully!", type = "message")
    }, error = function(e) {
      sc_log("ERROR", paste("Failed to load uploaded benchmark file:", e$message), "Upload")
      showNotification(paste("Upload error:", e$message), type = "error")
    })
  })
  
  
  output$ui_dataset_summary_download_btn <- renderUI({
    if (!is.null(rv$dataset_summary_df) && nrow(rv$dataset_summary_df) > 0) {
      downloadButton("download_dataset_summary_csv", "Export Properties (CSV)", class = "btn btn-sm btn-outline-secondary", icon = icon("file-csv"))
    }
  })
  
  filtered_dataset_summary_df <- reactive({
    req(rv$dataset_summary_df)
    df <- rv$dataset_summary_df
    filt <- input$opt_props_modality_filter %||% "all"
    if (filt == "rna") {
      df <- df[grepl("RNA", df$Modality, ignore.case = TRUE), , drop = FALSE]
    } else if (filt == "atac") {
      df <- df[grepl("ATAC", df$Modality, ignore.case = TRUE), , drop = FALSE]
    }
    df
  })

  output$ui_dataset_summary_kpis <- renderUI({
    if (is.null(rv$dataset_summary_df) || nrow(rv$dataset_summary_df) == 0) {
      return(p("Dataset properties will appear here once matrices are ingested.", style = "font-size: 0.85rem; color: #94A3B8; font-style: italic;"))
    }
    df <- filtered_dataset_summary_df()
    if (nrow(df) == 0) {
      return(div(class = "alert alert-warning py-2 px-3 my-2", style = "font-size: 0.84rem;",
                 icon("circle-info", class = "me-1"), "No datasets match the selected modality filter."))
    }
    filt <- input$opt_props_modality_filter %||% "all"
    ref_idx <- which(grepl("Reference", df$Role, ignore.case = TRUE))
    row_pick <- if (length(ref_idx) > 0) df[ref_idx[1], ] else df[1, ]
    
    feats_str <- if (filt == "all" && length(ref_idx) >= 2 && any(grepl("RNA", df$Modality, ignore.case = TRUE)) && any(grepl("ATAC", df$Modality, ignore.case = TRUE))) {
      rna_idx <- which(grepl("Reference", df$Role, ignore.case = TRUE) & grepl("RNA", df$Modality, ignore.case = TRUE))[1]
      atac_idx <- which(grepl("Reference", df$Role, ignore.case = TRUE) & grepl("ATAC", df$Modality, ignore.case = TRUE))[1]
      if (!is.na(rna_idx) && !is.na(atac_idx)) {
        paste0(df[rna_idx, "Features (P)"], " (RNA) + ", df[atac_idx, "Features (P)"], " (ATAC)")
      } else {
        row_pick[["Features (P)"]]
      }
    } else if (filt == "rna") {
      paste0(row_pick[["Features (P)"]], " (RNA)")
    } else if (filt == "atac") {
      paste0(row_pick[["Features (P)"]], " (ATAC)")
    } else {
      row_pick[["Features (P)"]]
    }
    
    fluidRow(
      column(2, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #1E3A8A; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", row_pick[["Cells (N)"]]),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Cells (N)"))),
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #0D9488; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", feats_str),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Features (P)"))),
      column(2, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #D97706; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", row_pick[["Sparsity"]]),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Sparsity (S%)"))),
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #6366F1; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", row_pick[["Cell Types (Groups)"]]),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Biological Groups (K)"))),
      column(2, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #0284C7; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", row_pick[["Batches"]]),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Batches (B)")))
    )
  })
  
  output$table_dataset_summary <- renderDT({
    req(rv$dataset_summary_df)
    df <- filtered_dataset_summary_df()
    bold_cols <- intersect(c("Dataset / Simulator", "Dataset", "Role", "Modality"), colnames(df))
    mono_cols <- intersect(c("Cells (N)", "Features (P)", "Sparsity", "Cell Types (Groups)", "Batches", "Median Lib Size", "Median Detected Features", "Mean Expression"), colnames(df))
    
    dt <- datatable(
      df,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = "ftip",
        autoWidth = TRUE
      ),
      rownames = FALSE,
      class = "compact stripe hover border"
    )
    if (length(bold_cols) > 0) {
      dt <- dt %>% formatStyle(columns = bold_cols, fontWeight = "bold")
    }
    if (length(mono_cols) > 0) {
      dt <- dt %>% formatStyle(columns = mono_cols, fontFamily = "monospace")
    }
    dt
  })
  output$download_dataset_summary_csv <- downloadHandler(
    filename = function() { paste0("dataset_properties_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv") },
    content = function(file) {
      req(rv$dataset_summary_df)
      utils::write.csv(rv$dataset_summary_df, file, row.names = FALSE)
    }
  )
  
  # ----------------------------------------------------------------------------
  # Data Hub: Main Panel UI & Dynamic Explorer
  # ----------------------------------------------------------------------------
  observeEvent(input$btn_jump_to_bubble, {
    nav_select("nav_active", "Comparative Bubble Matrix")
  })

  # Render Single-Cell Upload Workspace in Main Panel
  # ----------------------------------------------------------------------------
  # Data Hub: Standalone Ingestion & Property Inspection Loaders
  # ----------------------------------------------------------------------------
  load_and_inspect_unimodal_datasets <- function() {
    n_dsets <- input$num_uni_datasets %||% 1
    n_sims <- input$num_uni_sims %||% 1
    
    withProgress(message = "Reading Uploaded Datasets...", value = 0, {
      tryCatch({
        summary_rows <- list()
        all_mats <- list()
        refs_list <- list()
        sim_mats_list <- list()
        first_sim_mat <- NULL
        
        for (d in seq_len(n_dsets)) {
          ref_file <- input[[paste0("file_uni_ref_", d)]]
          dname_custom <- input[[paste0("uni_dset_name_", d)]]
          d_label <- if (!is.null(dname_custom) && trimws(dname_custom) != "") trimws(dname_custom) else sprintf("Dataset %d", d)
          
          ct_file <- input[[paste0("file_uni_celltypes_", d)]]
          bt_file <- input[[paste0("file_uni_batch_", d)]]
          ct_vec <- if (!is.null(ct_file) && nrow(ct_file) > 0) read_uploaded_labels(ct_file$datapath, ct_file$name) else NULL
          bt_vec <- if (!is.null(bt_file) && nrow(bt_file) > 0) read_uploaded_labels(bt_file$datapath, bt_file$name) else NULL
          
          if (!is.null(ct_vec)) {
            rv$cell_types_map[[d_label]] <- ct_vec
            rv$cell_types_map[[paste0(d_label, "_RNA")]] <- ct_vec
          }
          if (!is.null(bt_vec)) {
            rv$batches_map[[d_label]] <- bt_vec
            rv$batches_map[[paste0(d_label, "_RNA")]] <- bt_vec
          }
          
          if (!is.null(ref_file) && nrow(ref_file) > 0) {
            incProgress(0.5 / n_dsets, detail = sprintf("Reading %s reference...", d_label))
            ref_mat_d <- read_uploaded_matrix(ref_file$datapath, ref_file$name)
            if (inherits(ref_mat_d, "matrix") && !inherits(ref_mat_d, "dgCMatrix")) ref_mat_d <- Matrix::Matrix(ref_mat_d, sparse = TRUE)
            
            s_row <- extract_dataset_summary(
              ref_mat_d, role = "Biological Reference", method_name = paste0("Empirical Reference (", d_label, ")"),
              modality = "Single-Cell (Counts)", cell_types = ct_vec, batch_info = bt_vec
            )
            s_row$Dataset <- d_label
            summary_rows[[length(summary_rows) + 1]] <- s_row
            
            refs_list[[d_label]] <- ref_mat_d
            all_mats[[paste0("[", d_label, "] Empirical Reference")]] <- ref_mat_d
            
            if (d == 1) {
              rv$toy_ref <- ref_mat_d
              rv$cell_types <- ct_vec
              rv$batch <- bt_vec
            }
          }
          
          for (s in seq_len(n_sims)) {
            sim_name <- input[[paste0("uni_sim_name_", s)]]
            if (is.null(sim_name) || trimws(sim_name) == "") sim_name <- sprintf("Simulator %d", s)
            sim_file <- input[[paste0("file_uni_sim_", s, "_", d)]]
            
            if (!is.null(sim_file) && nrow(sim_file) > 0) {
              incProgress(0.5 / (n_dsets * n_sims), detail = sprintf("Reading %s on %s...", sim_name, d_label))
              sim_mat_d <- read_uploaded_matrix(sim_file$datapath, sim_file$name)
              if (inherits(sim_mat_d, "matrix") && !inherits(sim_mat_d, "dgCMatrix")) sim_mat_d <- Matrix::Matrix(sim_mat_d, sparse = TRUE)
              
              s_row <- extract_dataset_summary(
                sim_mat_d, role = "Simulated", method_name = paste0(sim_name, " (", d_label, ")"),
                modality = "Single-Cell (Counts)", cell_types = ct_vec, batch_info = bt_vec
              )
              s_row$Dataset <- d_label
              summary_rows[[length(summary_rows) + 1]] <- s_row
              
              all_mats[[paste0("[", d_label, "] Simulated: ", sim_name)]] <- sim_mat_d
              sim_mats_list[[paste0(sim_name, "_", d_label)]] <- sim_mat_d
              if (is.null(first_sim_mat)) first_sim_mat <- sim_mat_d
            }
          }
        }
        
        if (length(summary_rows) == 0) {
          showNotification("Please upload reference matrix or simulated matrix file(s) first.", type = "warning", duration = 4)
          return()
        }
        
        rv$dataset_summary_df <- do.call(rbind, summary_rows)
        rv$all_matrices <- all_mats
        rv$ref_matrices <- refs_list
        rv$sim_matrices <- sim_mats_list
        if (!is.null(first_sim_mat)) rv$toy_sim <- first_sim_mat
        
        showNotification(sprintf("Loaded properties for %d dataset/simulator matrices! Switched to Dataset properties.", length(summary_rows)), type = "message", duration = 5)
        nav_select("datahub_subtabs", "Dataset properties")
      }, error = function(e) {
        showNotification(paste("Error inspecting datasets:", e$message), type = "error", duration = 7)
      })
    })
  }

  load_and_inspect_multiomics_datasets <- function() {
    n_dsets <- input$num_multi_datasets %||% 1
    n_sims <- input$num_multi_sims %||% 1
    
    withProgress(message = "Reading Uploaded Multiomics Datasets...", value = 0, {
      tryCatch({
        summary_rows <- list()
        all_mats <- list()
        refs_list <- list()
        sim_mats_list <- list()
        first_sim_mat <- NULL
        
        for (d in seq_len(n_dsets)) {
          ref_rna_file <- input[[paste0("file_multi_ref_rna_", d)]]
          ref_atac_file <- input[[paste0("file_multi_ref_atac_", d)]]
          dname_custom <- input[[paste0("multi_dset_name_", d)]]
          d_label <- if (!is.null(dname_custom) && trimws(dname_custom) != "") trimws(dname_custom) else sprintf("Dataset %d", d)
          
          if (identical(input$opt_multi_pairing, "unpaired")) {
            ct_r_f <- input[[paste0("file_multi_celltypes_rna_", d)]] %||% input[[paste0("file_multi_celltypes_", d)]]
            ct_a_f <- input[[paste0("file_multi_celltypes_atac_", d)]]
            bt_r_f <- input[[paste0("file_multi_batch_rna_", d)]] %||% input[[paste0("file_multi_batch_", d)]]
            bt_a_f <- input[[paste0("file_multi_batch_atac_", d)]]
            ct_r_vec <- if (!is.null(ct_r_f) && nrow(ct_r_f) > 0) read_uploaded_labels(ct_r_f$datapath, ct_r_f$name) else NULL
            ct_a_vec <- if (!is.null(ct_a_f) && nrow(ct_a_f) > 0) read_uploaded_labels(ct_a_f$datapath, ct_a_f$name) else NULL
            bt_r_vec <- if (!is.null(bt_r_f) && nrow(bt_r_f) > 0) read_uploaded_labels(bt_r_f$datapath, bt_r_f$name) else NULL
            bt_a_vec <- if (!is.null(bt_a_f) && nrow(bt_a_f) > 0) read_uploaded_labels(bt_a_f$datapath, bt_a_f$name) else NULL
          } else {
            ct_file <- input[[paste0("file_multi_celltypes_", d)]]
            bt_file <- input[[paste0("file_multi_batch_", d)]]
            ct_r_vec <- if (!is.null(ct_file) && nrow(ct_file) > 0) read_uploaded_labels(ct_file$datapath, ct_file$name) else NULL
            ct_a_vec <- ct_r_vec
            bt_r_vec <- if (!is.null(bt_file) && nrow(bt_file) > 0) read_uploaded_labels(bt_file$datapath, bt_file$name) else NULL
            bt_a_vec <- bt_r_vec
          }
          
          if (!is.null(ct_r_vec)) {
            rv$cell_types_map[[paste0(d_label, "_RNA")]] <- ct_r_vec
            rv$cell_types_map[[d_label]] <- ct_r_vec
          }
          if (!is.null(ct_a_vec)) {
            rv$cell_types_map[[paste0(d_label, "_ATAC")]] <- ct_a_vec
          }
          if (!is.null(bt_r_vec)) {
            rv$batches_map[[paste0(d_label, "_RNA")]] <- bt_r_vec
            rv$batches_map[[d_label]] <- bt_r_vec
          }
          if (!is.null(bt_a_vec)) {
            rv$batches_map[[paste0(d_label, "_ATAC")]] <- bt_a_vec
          }
          if (d == 1) {
            rv$cell_types <- if (!is.null(ct_r_vec)) ct_r_vec else ct_a_vec
            rv$batch <- if (!is.null(bt_r_vec)) bt_r_vec else bt_a_vec
          }
          
          if (!is.null(ref_rna_file) && nrow(ref_rna_file) > 0) {
            incProgress(0.2 / n_dsets, detail = sprintf("Reading %s RNA reference...", d_label))
            ref_rna_d <- read_uploaded_matrix(ref_rna_file$datapath, ref_rna_file$name)
            if (inherits(ref_rna_d, "matrix") && !inherits(ref_rna_d, "dgCMatrix")) ref_rna_d <- Matrix::Matrix(ref_rna_d, sparse = TRUE)
            
            s_row <- extract_dataset_summary(
              ref_rna_d, role = "Biological Reference (RNA)", method_name = paste0("Empirical RNA (", d_label, ")"),
              modality = "RNA-seq", cell_types = ct_r_vec, batch_info = bt_r_vec
            )
            s_row$Dataset <- d_label
            summary_rows[[length(summary_rows) + 1]] <- s_row
            refs_list[[paste0(d_label, "_RNA")]] <- ref_rna_d
            all_mats[[paste0("[", d_label, "] Empirical RNA Reference")]] <- ref_rna_d
            if (d == 1) {
              rv$toy_ref <- ref_rna_d
            }
          }
          
          if (!is.null(ref_atac_file) && nrow(ref_atac_file) > 0) {
            incProgress(0.2 / n_dsets, detail = sprintf("Reading %s ATAC reference...", d_label))
            ref_atac_d <- read_uploaded_matrix(ref_atac_file$datapath, ref_atac_file$name)
            if (inherits(ref_atac_d, "matrix") && !inherits(ref_atac_d, "dgCMatrix")) ref_atac_d <- Matrix::Matrix(ref_atac_d, sparse = TRUE)
            
            s_row <- extract_dataset_summary(
              ref_atac_d, role = "Biological Reference (ATAC)", method_name = paste0("Empirical ATAC (", d_label, ")"),
              modality = "ATAC-seq", cell_types = ct_a_vec, batch_info = bt_a_vec
            )
            s_row$Dataset <- d_label
            summary_rows[[length(summary_rows) + 1]] <- s_row
            refs_list[[paste0(d_label, "_ATAC")]] <- ref_atac_d
            all_mats[[paste0("[", d_label, "] Empirical ATAC Reference")]] <- ref_atac_d
          }
          
          for (s in seq_len(n_sims)) {
            sim_name <- input[[paste0("multi_sim_name_", s)]]
            if (is.null(sim_name) || trimws(sim_name) == "") sim_name <- sprintf("Simulator %d", s)
            sim_rna_file  <- input[[paste0("file_multi_sim_rna_", s, "_", d)]]
            sim_atac_file <- input[[paste0("file_multi_sim_atac_", s, "_", d)]]
            
            if (!is.null(sim_rna_file) && nrow(sim_rna_file) > 0) {
              incProgress(0.3 / (n_dsets * n_sims), detail = sprintf("Reading %s RNA on %s...", sim_name, d_label))
              sim_rna_d <- read_uploaded_matrix(sim_rna_file$datapath, sim_rna_file$name)
              if (inherits(sim_rna_d, "matrix") && !inherits(sim_rna_d, "dgCMatrix")) sim_rna_d <- Matrix::Matrix(sim_rna_d, sparse = TRUE)
              
              s_row <- extract_dataset_summary(
                sim_rna_d, role = "Simulated (RNA)", method_name = paste0(sim_name, " RNA (", d_label, ")"),
                modality = "RNA-seq", cell_types = ct_r_vec, batch_info = bt_r_vec
              )
              s_row$Dataset <- d_label
              summary_rows[[length(summary_rows) + 1]] <- s_row
              all_mats[[paste0("[", d_label, "] Simulated RNA: ", sim_name)]] <- sim_rna_d
              sim_mats_list[[paste0(sim_name, "_RNA_", d_label)]] <- sim_rna_d
              if (is.null(first_sim_mat)) first_sim_mat <- sim_rna_d
            }
            
            if (!is.null(sim_atac_file) && nrow(sim_atac_file) > 0) {
              incProgress(0.3 / (n_dsets * n_sims), detail = sprintf("Reading %s ATAC on %s...", sim_name, d_label))
              sim_atac_d <- read_uploaded_matrix(sim_atac_file$datapath, sim_atac_file$name)
              if (inherits(sim_atac_d, "matrix") && !inherits(sim_atac_d, "dgCMatrix")) sim_atac_d <- Matrix::Matrix(sim_atac_d, sparse = TRUE)
              
              s_row <- extract_dataset_summary(
                sim_atac_d, role = "Simulated (ATAC)", method_name = paste0(sim_name, " ATAC (", d_label, ")"),
                modality = "ATAC-seq", cell_types = ct_a_vec, batch_info = bt_a_vec
              )
              s_row$Dataset <- d_label
              summary_rows[[length(summary_rows) + 1]] <- s_row
              all_mats[[paste0("[", d_label, "] Simulated ATAC: ", sim_name)]] <- sim_atac_d
              sim_mats_list[[paste0(sim_name, "_ATAC_", d_label)]] <- sim_atac_d
            }
          }
        }
        
        if (length(summary_rows) == 0) {
          showNotification("Please upload reference or simulated multiomics files first.", type = "warning", duration = 4)
          return()
        }
        
        rv$dataset_summary_df <- do.call(rbind, summary_rows)
        rv$all_matrices <- all_mats
        rv$ref_matrices <- refs_list
        rv$sim_matrices <- sim_mats_list
        if (!is.null(first_sim_mat)) rv$toy_sim <- first_sim_mat
        mode_data$multiomics$dataset_summary_df <- rv$dataset_summary_df
        mode_data$multiomics$all_matrices       <- all_mats
        mode_data$multiomics$ref_matrices       <- refs_list
        mode_data$multiomics$sim_matrices       <- sim_mats_list
        mode_data$multiomics$toy_ref            <- rv$toy_ref
        mode_data$multiomics$toy_sim            <- rv$toy_sim
        mode_data$multiomics$cell_types         <- rv$cell_types
        mode_data$multiomics$batch              <- rv$batch
        mode_data$multiomics$cell_types_map     <- rv$cell_types_map
        mode_data$multiomics$batches_map        <- rv$batches_map
        mode_data$multiomics$source_name        <- "Inspected Multiomics Datasets"
        sync_active_mode_data("multiomics")
        
        showNotification(sprintf("Loaded properties for %d multiomics matrices! Switched to Dataset properties.", length(summary_rows)), type = "message", duration = 5)
        nav_select("datahub_subtabs", "Dataset properties")
      }, error = function(e) {
        showNotification(paste("Error inspecting multiomics datasets:", e$message), type = "error", duration = 7)
      })
    })
  }

  observeEvent(input$btn_inspect_uni_data, { load_and_inspect_unimodal_datasets() })
  observeEvent(input$btn_inspect_uni_data_side, { load_and_inspect_unimodal_datasets() })
  observeEvent(input$btn_inspect_multi_data, { load_and_inspect_multiomics_datasets() })
  observeEvent(input$btn_inspect_multi_data_side, { load_and_inspect_multiomics_datasets() })
  observeEvent(input$btn_jump_to_upload_tab, { nav_select("datahub_subtabs", "Upload & Configure Data") })
  observeEvent(input$btn_jump_to_upload_tab2, { nav_select("datahub_subtabs", "Upload & Configure Data") })
  observeEvent(input$btn_back_to_upload, { nav_select("datahub_subtabs", "Upload & Configure Data") })
  observeEvent(input$btn_view_matrices_from_props, { nav_select("datahub_subtabs", "View Uploaded Data") })
  observeEvent(input$btn_goto_props_from_view, { nav_select("datahub_subtabs", "Dataset properties") })

  # Render Example Dataset Analysis Workspace in Main Panel
  render_demo_overview_workspace <- function() {
    has_demo_loaded <- !is.null(rv$benchmark_df)
    
    card(
      card_header(
        class = "bg-light",
        div(
          class = "d-flex justify-content-between align-items-center",
          tags$b("Example Dataset Analysis (Demo Datasets: 6 Simulators)", style = "color: #1E3A8A; font-size: 0.96rem;"),
          if (has_demo_loaded) {
            span(class = "badge bg-success text-white px-2 py-1", icon("check", class = "me-1"), "Loaded & Ready")
          } else {
            span(class = "badge bg-secondary text-white px-2 py-1", "Ready to Load")
          }
        )
      ),
      card_body(
        style = "padding: 24px;",
        div(
          class = "mb-3",
          p(
            "This section provides pre-computed demo datasets for exploring all scSimEval evaluation modules, diagnostic plots, and metric comparisons. ",
            "Because these demo datasets are pre-packaged for demonstration purposes, the methods are labeled as ",
            tags$strong("Simulator 1, Simulator 2, Simulator 3, Simulator 4, Simulator 5, and Simulator 6"),
            " evaluated against an empirical biological reference ground truth across ",
            tags$strong("62 canonical metrics"), " in 5 scientific categories.",
            style = "font-size: 0.88rem; color: #334155; line-height: 1.6;"
          )
        ),
        div(
          class = "row g-3 mb-4",
          div(
            class = "col-md-6",
            div(
              class = "p-3 rounded border h-100",
              style = "background: #F8FAFC;",
              h6(tags$b(icon("dna", class = "me-2 text-primary"), "2 Biological Reference Datasets"), style = "color: #1E293B; font-size: 0.88rem; margin-bottom: 6px;"),
              p("Empirical reference count matrices (Dataset 1 and Dataset 2) with cell-type and batch annotations serving as biological baseline.", style = "font-size: 0.80rem; color: #64748B; margin: 0;")
            )
          ),
          div(
            class = "col-md-6",
            div(
              class = "p-3 rounded border h-100",
              style = "background: #F8FAFC;",
              h6(tags$b(icon("layer-group", class = "me-2 text-primary"), "6 Demo Simulators"), style = "color: #1E293B; font-size: 0.88rem; margin-bottom: 6px;"),
              p("Simulated count matrices for each of the 2 reference datasets (Simulator 1 to Simulator 6) evaluated across 62 canonical metrics.", style = "font-size: 0.80rem; color: #64748B; margin: 0;")
            )
          )
        ),
        if (!has_demo_loaded) {
          div(
            class = "text-center py-2",
            actionButton("btn_empty_load_demo", "Load Example Benchmark (2 Datasets x 6 Simulators)", class = "btn btn-primary px-4 py-2", icon = icon("play"), style = "font-weight: 600; font-size: 0.95rem;")
          )
        } else {
          div(
            class = "alert alert-success d-flex justify-content-between align-items-center flex-wrap gap-2 p-3 mb-0",
            style = "background: #F0FDF4; border: 1px solid #BBF7D0; border-radius: 6px;",
            div(
              icon("check-circle", class = "text-success me-2", style = "font-size: 1.25rem;"),
              span(style = "color: #166534; font-weight: 600; font-size: 0.90rem;",
                   "Demo benchmark loaded: 2 reference datasets (Dataset 1, Dataset 2) and 6 simulators (Simulator 1 to Simulator 6 across 62 metrics)."),
              div(style = "color: #15803D; font-size: 0.82rem; margin-top: 2px;", 
                  "Inspect matrices in 'View Uploaded Data' or examine properties and rankings in 'Dataset properties':")
            ),
            div(
              class = "d-flex gap-2 flex-wrap",
              actionButton("btn_demo_goto_view", "View Uploaded Data", class = "btn btn-sm btn-outline-success", icon = icon("magnifying-glass")),
              actionButton("btn_demo_goto_props", "Dataset properties", class = "btn btn-sm btn-outline-success", icon = icon("table-cells")),
              actionButton("btn_demo_goto_bubble", "Comparative Bubble Matrix →", class = "btn btn-sm btn-success", icon = icon("chart-simple"))
            )
          )
        }
      )
    )
  }

  # Render Single-Cell Upload Workspace in Main Panel
  render_unimodal_upload_workspace <- function() {
    n_dsets <- input$num_uni_datasets %||% 1
    n_sims <- input$num_uni_sims %||% 1
    
    div(
      div(
        class = "mb-3 p-3 rounded",
        style = "background: linear-gradient(135deg, #EFF6FF 0%, #DBEAFE 100%); border: 1px solid #BFDBFE;",
        div(
          class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
          div(
            h5(tags$b("Single-Cell Simulation Benchmarking Workspace"), style = "color: #1E3A8A; margin-bottom: 4px; font-size: 1.12rem;"),
            p("Upload biological reference and simulated count matrices (scRNA-seq or scATAC-seq) along with computational benchmarks (Runtime, Peak RAM).",
              style = "color: #1E40AF; margin: 0; font-size: 0.84rem;")
          ),
          div(
            class = "d-flex gap-2 align-items-center flex-wrap",
            span(class = "badge bg-primary text-white px-3 py-2", style = "font-size: 0.85rem;",
                 sprintf("%d Reference Dataset%s | %d Simulator%s", n_dsets, if (n_dsets > 1) "s" else "", n_sims, if (n_sims > 1) "s" else "")),
            uiOutput("ui_uni_cores_badge", inline = TRUE)
          )
        )
      ),
      
            # Metric Capacity Callout Banner
      div(
        class = "alert alert-info py-2 px-3 mb-3 border-info-subtle",
        style = "background: #F0F9FF; border-left: 4px solid #0284C7; font-size: 0.83rem; color: #0369A1;",
        div(class = "d-flex align-items-center justify-content-between flex-wrap gap-2",
          div(
            style = "color: #0369A1; line-height: 1.55;",
            tags$b(icon("circle-info", class = "me-1 text-info"), "Exact Single-Cell Metric Measurement Capacity: ", style = "color: #0C4A6E;"),
            "Count matrix alone measures ", tags$b("24 metrics (scRNA)", style = "color: #0C4A6E;"), " / ", tags$b("22 metrics (scATAC)", style = "color: #0C4A6E;"), " (Distributional properties, correlations, zero-inflation curves, trajectory). ",
            "Uploading ", tags$span(class = "badge bg-primary-subtle text-primary border border-primary-subtle", "+ Cell Types"), 
            " unlocks ", tags$b("44 metrics (scRNA)", style = "color: #0C4A6E;"), " / ", tags$b("36 metrics (scATAC)", style = "color: #0C4A6E;"), " (clustering concordance, ARI, NMI, SimBench biological signal fidelity). ",
            "Uploading ", tags$span(class = "badge bg-warning-subtle text-warning border border-warning-subtle", "+ Batch Info"), 
            " unlocks ", tags$b("51 metrics (scRNA)", style = "color: #0C4A6E;"), " / ", tags$b("43 metrics (scATAC)", style = "color: #0C4A6E;"), " (batch mixing, kBET, LISI, CMS, Shannon entropy)."
          ),
          actionLink(
            "btn_show_metadata_format_hint",
            span(icon("circle-question", class = "me-1"), "View Formats & Full Breakdown →"),
            class = "btn btn-sm btn-outline-primary fw-bold shadow-sm",
            style = "font-size: 0.78rem; white-space: nowrap; border-color: #0284C7; color: #0284C7 !important; background: #FFFFFF;"
          )
        )
      ),
      # Card 1: Common Biological Reference Single-Cell Datasets (1 to N)
      card(
        class = "mb-3",
        card_header(
          class = "bg-light",
          div(
            class = "d-flex justify-content-between align-items-center",
            tags$b("1. Original Reference Datasets (Common ground truth across all simulators)", style = "color: #1E3A8A; font-size: 0.96rem;"),
            span(class = "badge bg-info text-white", "Uploaded Once")
          )
        ),
        card_body(
          p("Upload empirical reference count matrices and optional annotations for each dataset. All simulators below will be benchmarked against these references:",
            style = "font-size: 0.82rem; color: #64748B; margin-bottom: 12px;"),
          uiOutput("ui_uni_ref_cards")
        )
      ),
      
      # Card 2: Simulators & Simulated Single-Cell Datasets (1 to M)
      card(
        class = "mb-3",
        card_header(
          class = "bg-light",
          div(
            class = "d-flex justify-content-between align-items-center",
            tags$b("2. Simulators & Simulated Datasets (Click to expand/hide individual simulators)", style = "color: #1E3A8A; font-size: 0.96rem;"),
            span(class = "badge bg-primary", sprintf("%d Simulators", n_sims))
          )
        ),
        card_body(
          p("For each simulator, upload the simulated count matrix corresponding to each dataset above, along with observed runtime and Peak RAM (MiB):",
            style = "font-size: 0.82rem; color: #64748B; margin-bottom: 12px;"),
          uiOutput("ui_uni_sim_cards")
        )
      ),
      
      # Scientific Callout Banner
      div(
        class = "alert alert-secondary py-2 px-3 mt-3 mb-2 rounded",
        style = "font-size: 0.80rem; border: 1px dashed #CBD5E1;",
        icon("file-code", class = "me-2 text-primary"),
        tags$strong("Reproducible Analysis:"),
        " Click ", tags$b("Download R Script"),
        " to export an executable R script pre-configured with all above dataset names, paths, and benchmark parameters."
      ),
      
      # Bottom Action Bar
      div(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2 mt-2 p-3 bg-light border rounded",
        div(
          class = "d-flex align-items-center gap-3 flex-wrap",
          checkboxInput("chk_append_uni", "Add to current benchmark results", value = FALSE),
          uiOutput("ui_uni_cores_badge_bottom", inline = TRUE)
        ),
        div(
          class = "d-flex gap-2 flex-wrap",
          downloadButton("dl_hpc_r_script_uni", "Download R Script", class = "btn btn-outline-primary", icon = icon("file-code")),
          actionButton("btn_inspect_uni_data", "Inspect & View Dataset Properties", class = "btn btn-outline-secondary px-3 py-2", icon = icon("table-cells")),
          actionButton("btn_run_uni_eval_main", "Calculate Single-Cell Metrics", class = "btn btn-primary px-4 py-2", icon = icon("play"), style = "font-weight: 700; font-size: 0.95rem;")
        )
      )
    )
  }

  # Render Multiomics Upload Workspace in Main Panel
  render_multiomics_upload_workspace <- function() {
    n_dsets <- input$num_multi_datasets %||% 1
    n_sims <- input$num_multi_sims %||% 1
    
    div(
      div(
        class = "mb-3 p-3 rounded",
        style = "background: linear-gradient(135deg, #EEF2FF 0%, #E0E7FF 100%); border: 1px solid #C7D2FE;",
        div(
          class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
          div(
            h5(tags$b("Multiomics Simulation Benchmarking Workspace"), style = "color: #312E81; margin-bottom: 4px; font-size: 1.12rem;"),
            p("Upload paired reference RNA & ATAC datasets once as ground truth. Then, under each simulator, upload simulated RNA & ATAC matrices.",
              style = "color: #4338CA; margin: 0; font-size: 0.84rem;")
          ),
          div(
            class = "d-flex gap-2 align-items-center flex-wrap",
            span(class = "badge text-white px-3 py-2", style = "background: #4F46E5; font-size: 0.85rem;",
                 sprintf("%d Reference Dataset%s | %d Simulator%s", n_dsets, if (n_dsets > 1) "s" else "", n_sims, if (n_sims > 1) "s" else "")),
            uiOutput("ui_multi_cores_badge", inline = TRUE)
          )
        )
      ),
      
            # Metric Capacity Callout Banner
      div(
        class = "alert alert-primary py-2 px-3 mb-3 border-indigo-subtle",
        style = "background: #EEF2FF; border-left: 4px solid #4F46E5; font-size: 0.83rem; color: #3730A3;",
        div(class = "d-flex align-items-center justify-content-between flex-wrap gap-2",
          div(
            style = "color: #3730A3; line-height: 1.55;",
            tags$b(icon("circle-info", class = "me-1 text-indigo"), "Exact Multiomics Metric Measurement Capacity: ", style = "color: #1E1B4B;"),
            "Paired RNA + ATAC matrices alone evaluate ", tags$b("36 / 62 metrics", style = "color: #1E1B4B;"), " (unimodal properties for both layers + Category VII Cross-Modal Coupling, FOSCTTM alignment, Match@1, Match@5, and peak-to-gene linkages). ",
            "Including ", tags$span(class = "badge bg-primary-subtle text-primary border border-primary-subtle", "+ Cell Types"), " unlocks ",
            tags$b("55 / 62 metrics", style = "color: #1E1B4B;"), ". Adding ",
            tags$span(class = "badge bg-warning-subtle text-warning border border-warning-subtle", "+ Batch Info"), 
            " unlocks ", tags$b("all 62 / 62 metrics across all 8 categories", style = "color: #1E1B4B;"), ". ",
            "In Unpaired mode: 30 / 62 metrics (counts alone) and up to 52 / 62 metrics (with cell types and batches)."
          ),
          actionLink(
            "btn_show_metadata_format_multi_hint",
            span(icon("circle-question", class = "me-1"), "View Formats & Full Breakdown →"),
            class = "btn btn-sm btn-outline-primary fw-bold shadow-sm",
            style = "font-size: 0.78rem; white-space: nowrap; border-color: #4F46E5; color: #4F46E5 !important; background: #FFFFFF;"
          )
        )
      ),
      # Card 1: Common Biological Reference Multiomics Datasets (1 to N)
      card(
        class = "mb-3",
        card_header(
          class = "bg-light",
          div(
            class = "d-flex justify-content-between align-items-center",
            tags$b("1. Original Reference Multiomics Datasets (Common ground truth across all simulators)", style = "color: #312E81; font-size: 0.96rem;"),
            span(class = "badge bg-info text-white", "Uploaded Once")
          )
        ),
        card_body(
          p("Upload empirical reference RNA and ATAC count matrices and optional cell metadata for each dataset:",
            style = "font-size: 0.82rem; color: #64748B; margin-bottom: 12px;"),
          uiOutput("ui_multi_ref_cards")
        )
      ),
      
      # Card 2: Simulators & Simulated Multiomics Datasets (1 to M)
      card(
        class = "mb-3",
        card_header(
          class = "bg-light",
          div(
            class = "d-flex justify-content-between align-items-center",
            tags$b("2. Multiomics Simulators & Simulated Outputs (Click to expand/hide individual simulators)", style = "color: #312E81; font-size: 0.96rem;"),
            span(class = "badge text-white", style = "background: #4F46E5;", sprintf("%d Simulators", n_sims))
          )
        ),
        card_body(
          p("For each simulator, upload simulated RNA and ATAC count matrices corresponding to each dataset above, along with runtime and Peak RAM (MiB):",
            style = "font-size: 0.82rem; color: #64748B; margin-bottom: 12px;"),
          uiOutput("ui_multi_sim_cards")
        )
      ),
      
      # Scientific Callout Banner
      div(
        class = "alert alert-secondary py-2 px-3 mt-3 mb-2 rounded",
        style = "font-size: 0.80rem; border: 1px dashed #CBD5E1;",
        icon("file-code", class = "me-2 text-primary"),
        tags$strong("Reproducible Analysis:"),
        " Click ", tags$b("Download R Script"),
        " to export an executable R script pre-configured with all above multiomics dataset names, paths, and benchmark parameters."
      ),
      
      # Bottom Action Bar
      div(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2 mt-2 p-3 bg-light border rounded",
        div(
          class = "d-flex align-items-center gap-3 flex-wrap",
          checkboxInput("chk_append_multi", "Add to current benchmark results", value = FALSE),
          uiOutput("ui_multi_cores_badge_bottom", inline = TRUE)
        ),
        div(
          class = "d-flex gap-2 flex-wrap",
          downloadButton("dl_hpc_r_script_multi", "Download R Script", class = "btn btn-outline-primary", icon = icon("file-code")),
          actionButton("btn_inspect_multi_data", "Inspect & View Dataset Properties", class = "btn btn-outline-secondary px-3 py-2", icon = icon("table-cells")),
          actionButton("btn_run_multi_eval_main", "Calculate Multiomics Metrics", class = "btn btn-primary px-4 py-2", icon = icon("dna"), style = "font-weight: 700; font-size: 0.95rem;")
        )
      )
    )
  }

  # Render Data Hub Main Panel with 3 Unified Subtabs
  output$ui_datahub_main_panel <- renderUI({
    mode <- input$opt_data_mode %||% "demo"
    
    navset_card_tab(
      id = "datahub_subtabs",
      selected = "Upload & Configure Data",
      
      # Subtab 1: Upload & Configure Data
      nav_panel(
        title = span(icon("sliders", class = "me-1"), "Upload & Configure Data"),
        value = "Upload & Configure Data",
        div(
          style = "padding: 8px 4px;",
          if (mode == "demo") {
            render_demo_overview_workspace()
          } else if (mode == "multiomics") {
            render_multiomics_upload_workspace()
          } else if (mode == "upload_bench") {
            card(
              card_header(tags$b("Upload Saved Benchmark Results (.rds / .csv)", style = "color: #1E3A8A; font-size: 0.95rem;")),
              card_body(
                style = "padding: 30px;",
                p("Load a previously saved benchmark result file (.rds or .csv) from scSimEval to inspect visualizations and explore metrics.", style = "color: #64748B;"),
                fileInput("file_bench_upload_main", "Select Saved File(s) (.rds / .csv):", multiple = TRUE, accept = c(".rds", ".csv")),
                actionButton("btn_load_uploaded_bench_main", "Load Saved Results", class = "btn btn-primary px-4 py-2", icon = icon("folder-open"))
              )
            )
          } else {
            render_unimodal_upload_workspace()
          }
        )
      ),
      
      # Subtab 2: View Uploaded Data
      nav_panel(
        title = span(icon("magnifying-glass", class = "me-1"), "View Uploaded Data"),
        value = "View Uploaded Data",
        div(
          style = "padding: 8px 4px;",
          uiOutput("ui_datahub_view_data_subtab")
        )
      ),
      
      # Subtab 3: Dataset properties
      nav_panel(
        title = span(icon("table-cells", class = "me-1"), "Dataset properties"),
        value = "Dataset properties",
        div(
          style = "padding: 8px 4px;",
          uiOutput("ui_datahub_properties_subtab")
        )
      )
    )
  })

  # Subtab 3 Renderer: Dataset properties
  output$ui_datahub_properties_subtab <- renderUI({
    mode <- input$opt_data_mode %||% "demo"
    if (is.null(rv$dataset_summary_df) || nrow(rv$dataset_summary_df) == 0) {
      return(
        div(
          class = "text-center p-5 bg-light border rounded",
          style = "max-width: 650px; margin: 30px auto;",
          icon("table-cells", class = "text-muted mb-3", style = "font-size: 2.8rem; opacity: 0.6;"),
          h5("No Dataset Properties Available Yet", style = "font-weight: 700; color: #1E293B; margin-bottom: 8px;"),
          if (mode == "demo") {
            tagList(
              p("Load the pre-computed demo datasets (2 reference datasets and 6 simulators) to view comprehensive cell counts, feature dimensions, sparsity, and library depths across simulated outputs and empirical references.",
                style = "font-size: 0.90rem; color: #64748B; margin-bottom: 20px; line-height: 1.5;"),
              actionButton("btn_empty_load_demo_props", "Load Example Benchmark (2 Datasets x 6 Simulators)", class = "btn btn-primary px-4 py-2", icon = icon("play"))
            )
          } else if (mode == "multiomics") {
            tagList(
              p("Upload your paired or unpaired reference and simulated RNA & ATAC datasets in the ",
                tags$strong("Upload & Configure Data"),
                " subtab, then click ", tags$b("Inspect & View Dataset Properties"),
                " or ", tags$b("Calculate Multiomics Metrics"), " to view comprehensive cell counts, feature dimensions, sparsity, and library depths across all datasets.",
                style = "font-size: 0.90rem; color: #64748B; margin-bottom: 20px; line-height: 1.5;"),
              actionButton("btn_jump_to_upload_tab", "Go to Upload & Configure Data", class = "btn btn-primary px-4 py-2", icon = icon("sliders"))
            )
          } else if (mode == "upload_bench") {
            tagList(
              p("Upload a previously saved benchmark result file (.rds or .csv) in the ",
                tags$strong("Upload & Configure Data"),
                " subtab to inspect results and dataset properties.",
                style = "font-size: 0.90rem; color: #64748B; margin-bottom: 20px; line-height: 1.5;"),
              actionButton("btn_jump_to_upload_tab", "Go to Upload & Configure Data", class = "btn btn-primary px-4 py-2", icon = icon("sliders"))
            )
          } else {
            tagList(
              p("Upload your biological reference and simulated single-cell datasets in the ",
                tags$strong("Upload & Configure Data"),
                " subtab, then click ", tags$b("Inspect & View Dataset Properties"),
                " or ", tags$b("Calculate Single-Cell Metrics"), " to view comprehensive cell counts, feature dimensions, sparsity, and library depths across all datasets.",
                style = "font-size: 0.90rem; color: #64748B; margin-bottom: 20px; line-height: 1.5;"),
              actionButton("btn_jump_to_upload_tab", "Go to Upload & Configure Data", class = "btn btn-primary px-4 py-2", icon = icon("sliders"))
            )
          }
        )
      )
    }
    
    tagList(
      uiOutput("ui_status_banner"),
      div(
        class = "dataset-summary-box mb-3",
        div(
          class = "d-flex justify-content-between align-items-center flex-wrap gap-2 mt-2 mb-2",
          div(
            h6(tags$b("Dataset Properties Summary"), style = "color: #0F172A; margin: 0; font-size: 0.98rem;"),
            p("Number of cells, features, sparsity, and sequencing depth across empirical reference and simulated datasets:",
              style = "font-size: 0.82rem; color: #64748B; margin: 0;")
          ),
          div(
            class = "d-flex gap-2 flex-wrap",
            downloadButton("download_dataset_summary_csv", "Export Table (CSV)", class = "btn btn-sm btn-outline-primary", icon = icon("file-csv")),
            actionButton("btn_view_matrices_from_props", "View Uploaded Data", class = "btn btn-sm btn-outline-secondary", icon = icon("magnifying-glass")),
            if (mode != "demo") actionButton("btn_back_to_upload", "Upload & Configure Data", class = "btn btn-sm btn-outline-secondary", icon = icon("sliders"))
          )
        ),
        if (mode == "multiomics" || (!is.null(rv$dataset_summary_df) && any(grepl("RNA", rv$dataset_summary_df$Modality, ignore.case = TRUE)) && any(grepl("ATAC", rv$dataset_summary_df$Modality, ignore.case = TRUE)))) {
          div(
            style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 8px 14px; margin-top: 8px; margin-bottom: 12px; display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 10px;",
            div(
              style = "font-weight: 700; color: #1E293B; font-size: 0.86rem; display: flex; align-items: center; gap: 6px;",
              icon("filter", class = "text-primary"),
              "Filter Modality:"
            ),
            radioButtons(
              "opt_props_modality_filter", NULL,
              choices = c("All Modalities (scRNA + scATAC)" = "all", "scRNA-seq" = "rna", "scATAC-seq" = "atac"),
              selected = "all", inline = TRUE
            )
          )
        },
        uiOutput("ui_dataset_summary_kpis"),
        div(style = "margin-top: 14px; border: 1px solid #E2E8F0; border-radius: 6px; padding: 6px; background: #FFFFFF;",
            DTOutput("table_dataset_summary"))
      ),
      uiOutput("ui_datahub_annotations_breakdown")
    )
  })

  # Subtab 2 Renderer: View Uploaded Data
  output$ui_datahub_view_data_subtab <- renderUI({
    has_mats <- (!is.null(rv$all_matrices) && length(rv$all_matrices) > 0) || !is.null(rv$toy_ref) || (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0)
    mode <- input$opt_data_mode %||% "demo"
    if (!has_mats) {
      return(
        div(
          class = "text-center p-5 bg-light border rounded",
          style = "max-width: 650px; margin: 30px auto;",
          icon("magnifying-glass", class = "text-muted mb-3", style = "font-size: 2.8rem; opacity: 0.6;"),
          h5("No Matrices Available for Inspection", style = "font-weight: 700; color: #1E293B; margin-bottom: 8px;"),
          if (mode == "demo") {
            tagList(
              p("Load the pre-computed demo datasets (2 reference datasets and 6 simulators) to inspect count matrices (Head values, cell statistics, and feature statistics) across empirical references and simulated datasets.",
                style = "font-size: 0.90rem; color: #64748B; margin-bottom: 20px; line-height: 1.5;"),
              actionButton("btn_empty_load_demo_view", "Load Example Benchmark (2 Datasets x 6 Simulators)", class = "btn btn-primary px-4 py-2", icon = icon("play"))
            )
          } else if (mode == "multiomics") {
            tagList(
              p("Upload your reference and simulated RNA & ATAC count matrices in the ", tags$strong("Upload & Configure Data"),
                " subtab, then click ", tags$b("Inspect & View Dataset Properties"),
                " to inspect matrix heads, library size distributions, and gene/peak statistics.",
                style = "font-size: 0.90rem; color: #64748B; margin-bottom: 20px; line-height: 1.5;"),
              actionButton("btn_jump_to_upload_tab2", "Go to Upload & Configure Data", class = "btn btn-primary px-4 py-2", icon = icon("sliders"))
            )
          } else if (mode == "upload_bench") {
            tagList(
              p("Upload a saved benchmark result file in the ", tags$strong("Upload & Configure Data"),
                " subtab to explore evaluated metrics and visualizations.",
                style = "font-size: 0.90rem; color: #64748B; margin-bottom: 20px; line-height: 1.5;"),
              actionButton("btn_jump_to_upload_tab2", "Go to Upload & Configure Data", class = "btn btn-primary px-4 py-2", icon = icon("sliders"))
            )
          } else {
            tagList(
              p("Upload count matrices in the ", tags$strong("Upload & Configure Data"),
                " subtab, then click ", tags$b("Inspect & View Dataset Properties"),
                " to inspect raw count head values, library size distributions, and gene detection statistics across your single-cell datasets.",
                style = "font-size: 0.90rem; color: #64748B; margin-bottom: 20px; line-height: 1.5;"),
              actionButton("btn_jump_to_upload_tab2", "Go to Upload & Configure Data", class = "btn btn-primary px-4 py-2", icon = icon("sliders"))
            )
          }
        )
      )
    }
    is_multi_view <- identical(mode, "multiomics") || (has_mats && any(grepl("RNA", names(rv$all_matrices %||% rv$sim_matrices), ignore.case = TRUE)) && any(grepl("ATAC", names(rv$all_matrices %||% rv$sim_matrices), ignore.case = TRUE)))
    
    div(
      div(
        style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 12px 16px; margin-bottom: 14px;",
        fluidRow(
          if (is_multi_view) {
            column(
              3,
              radioButtons(
                "opt_inspect_modality_filter", tags$b(icon("filter", class = "me-1 text-primary"), "Modality:"),
                choices = c("All" = "all", "scRNA-seq" = "rna", "scATAC-seq" = "atac"),
                inline = TRUE, selected = "all"
              )
            )
          },
          column(
            if (is_multi_view) 4 else 5,
            uiOutput("ui_inspect_dataset_selector")
          ),
          column(
            if (is_multi_view) 3 else 4,
            radioButtons(
              "opt_inspect_view_mode", "View Format:",
              choices = c("Raw Values (Head)" = "matrix_head", "Cell Statistics" = "cell_stats", "Feature Statistics" = "feature_stats"),
              inline = TRUE, selected = "matrix_head"
            )
          ),
          column(
            2,
            div(class = "text-end mt-3",
                actionButton("btn_goto_props_from_view", "Dataset properties", class = "btn btn-sm btn-outline-secondary me-1", icon = icon("table-cells")),
                downloadButton("download_inspect_data_csv", "Export (CSV)", class = "btn btn-sm btn-outline-primary", icon = icon("download"))
            )
          )
        )
      ),
      uiOutput("ui_inspect_matrix_kpis"),
      div(style = "margin-top: 10px; border: 1px solid #E2E8F0; border-radius: 6px; padding: 6px; background: #FFFFFF;",
          DTOutput("table_matrix_inspect"))
    )
  })
  output$ui_datahub_annotations_breakdown <- renderUI({
    has_ct <- !is.null(rv$cell_types) && length(rv$cell_types) > 0
    has_bt <- !is.null(rv$batch) && length(rv$batch) > 0
    if (!has_ct && !has_bt) return(NULL)
    
    fluidRow(
      style = "margin-top: 14px;",
      if (has_ct) {
        ct_tbl <- as.data.frame(table(Cell_Type = rv$cell_types))
        ct_tbl$Percentage <- sprintf("%.1f%%", ct_tbl$Freq / sum(ct_tbl$Freq) * 100)
        colnames(ct_tbl) <- c("Cell Type / Cluster", "Cell Count", "Percentage")
        column(
          if (has_bt) 6 else 12,
          div(
            style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 12px;",
            h6(tags$b(icon("layer-group", class = "me-1"), " Biological Group Distribution"), style = "color: #1E3A8A; margin-bottom: 8px; font-size: 0.88rem;"),
            renderTable(ct_tbl, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "s")
          )
        )
      },
      if (has_bt) {
        bt_tbl <- as.data.frame(table(Batch = rv$batch))
        bt_tbl$Percentage <- sprintf("%.1f%%", bt_tbl$Freq / sum(bt_tbl$Freq) * 100)
        colnames(bt_tbl) <- c("Technical Batch", "Cell Count", "Percentage")
        column(
          if (has_ct) 6 else 12,
          div(
            style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 12px;",
            h6(tags$b(icon("boxes-stacked", class = "me-1"), " Technical Batch Distribution"), style = "color: #0D9488; margin-bottom: 8px; font-size: 0.88rem;"),
            renderTable(bt_tbl, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "s")
          )
        )
      }
    )
  })

  output$ui_inspect_dataset_selector <- renderUI({
    choices <- list()
    filt <- input$opt_inspect_modality_filter %||% "all"
    
    if (!is.null(rv$all_matrices) && length(rv$all_matrices) > 0) {
      for (mname in names(rv$all_matrices)) {
        if (filt == "rna" && !grepl("RNA", mname, ignore.case = TRUE)) next
        if (filt == "atac" && !grepl("ATAC", mname, ignore.case = TRUE)) next
        choices[[mname]] <- mname
      }
    } else {
      if (!is.null(rv$ref_matrices) && length(rv$ref_matrices) > 0) {
        for (rname in names(rv$ref_matrices)) {
          if (filt == "rna" && !grepl("RNA", rname, ignore.case = TRUE)) next
          if (filt == "atac" && !grepl("ATAC", rname, ignore.case = TRUE)) next
          choices[[paste0("Reference: ", rname)]] <- paste0("ref_", rname)
        }
      } else if (!is.null(rv$toy_ref)) {
        if (filt != "atac") choices[["Empirical Reference (Real Cells)"]] <- "ref"
      }
      if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) {
        for (sname in names(rv$sim_matrices)) {
          if (filt == "rna" && !grepl("RNA", sname, ignore.case = TRUE)) next
          if (filt == "atac" && !grepl("ATAC", sname, ignore.case = TRUE)) next
          choices[[paste0("Simulated: ", sname)]] <- paste0("sim_", sname)
        }
      } else if (!is.null(rv$toy_sim)) {
        if (filt != "atac") choices[["Simulated Matrix"]] <- "toy_sim"
      }
    }
    if (length(choices) == 0) choices[["No Matching Matrix Available"]] <- "none"
    selectInput("sel_inspect_matrix_choice", "Select Uploaded Dataset / Matrix:", choices = choices, selected = choices[[1]])
  })

  active_inspect_data <- reactive({
    choice <- input$sel_inspect_matrix_choice
    if (is.null(choice) || identical(choice, "none")) return(NULL)
    mode <- input$opt_inspect_view_mode %||% "matrix_head"
    
    mat <- NULL
    if (!is.null(choice) && !is.null(rv$all_matrices) && choice %in% names(rv$all_matrices)) {
      mat <- rv$all_matrices[[choice]]
    } else if (identical(choice, "ref") && !is.null(rv$toy_ref)) {
      mat <- rv$toy_ref
    } else if (identical(choice, "toy_sim") && !is.null(rv$toy_sim)) {
      mat <- rv$toy_sim
    } else if (!is.null(choice) && startsWith(choice, "sim_") && !is.null(rv$sim_matrices)) {
      sname <- sub("^sim_", "", choice)
      mat <- rv$sim_matrices[[sname]]
    } else if (!is.null(choice) && !is.null(rv$sim_matrices) && choice %in% names(rv$sim_matrices)) {
      mat <- rv$sim_matrices[[choice]]
    } else if (!is.null(rv$toy_ref)) {
      mat <- rv$toy_ref
    } else if (!is.null(rv$toy_sim)) {
      mat <- rv$toy_sim
    }
    if (is.null(mat)) return(NULL)
    
    if (mode %in% c("matrix_head", "raw", "lognorm")) {
      nr <- min(50, nrow(mat))
      nc <- min(25, ncol(mat))
      sub_mat <- as.matrix(mat[1:nr, 1:nc, drop = FALSE])
      if (mode == "lognorm") sub_mat <- round(log1p(sub_mat), 3)
      df <- as.data.frame(sub_mat)
      feat_names <- rownames(sub_mat)
      if (is.null(feat_names)) feat_names <- paste0("Feature_", seq_len(nr))
      df <- cbind(Feature_ID = feat_names, df)
      rownames(df) <- NULL
      list(df = df, mat = mat, type = "matrix_head")
    } else if (mode %in% c("feature_stats", "stats")) {
      mean_expr <- Matrix::rowMeans(mat)
      det_rate <- Matrix::rowMeans(mat > 0)
      if (inherits(mat, "dgCMatrix")) {
        sq_means <- Matrix::rowMeans(mat^2)
        vars <- (sq_means - mean_expr^2) * (ncol(mat) / max(1, ncol(mat) - 1))
        vars[vars < 0] <- 0
      } else {
        vars <- apply(mat, 1, stats::var)
      }
      feat_names <- rownames(mat)
      if (is.null(feat_names)) feat_names <- paste0("Feature_", seq_len(nrow(mat)))
      df <- data.frame(
        Feature = feat_names,
        Mean_Expression = round(as.numeric(mean_expr), 4),
        Variance = round(as.numeric(vars), 4),
        Detection_Rate_Pct = round(as.numeric(det_rate) * 100, 2),
        Dropout_Rate_Pct = round((1 - as.numeric(det_rate)) * 100, 2),
        stringsAsFactors = FALSE
      )
      df <- df[order(-df$Mean_Expression), ]
      rownames(df) <- NULL
      list(df = df, mat = mat, type = "feature_stats")
    } else {
      lib_sizes <- Matrix::colSums(mat)
      det_genes <- Matrix::colSums(mat > 0)
      cell_names <- colnames(mat)
      if (is.null(cell_names)) cell_names <- paste0("Cell_", seq_len(ncol(mat)))
      df <- data.frame(
        Cell_ID = cell_names,
        Library_Size = as.numeric(lib_sizes),
        Detected_Features = as.integer(det_genes),
        Sparsity_Pct = round((1 - as.numeric(det_genes) / max(1, nrow(mat))) * 100, 2),
        stringsAsFactors = FALSE
      )
      if (!is.null(rv$cell_types) && length(rv$cell_types) == ncol(mat)) {
        df$Cell_Type <- rv$cell_types
      }
      if (!is.null(rv$batch) && length(rv$batch) == ncol(mat)) {
        df$Batch <- rv$batch
      }
      rownames(df) <- NULL
      list(df = df, mat = mat, type = "cell_stats")
    }
  })
  output$ui_inspect_matrix_kpis <- renderUI({
    res <- active_inspect_data()
    if (is.null(res) || is.null(res$mat)) return(NULL)
    mat <- res$mat
    n_cells <- ncol(mat)
    n_feats <- nrow(mat)
    sparsity_val <- round(mean(mat == 0, na.rm = TRUE) * 100, 1)
    med_depth <- round(stats::median(Matrix::colSums(mat)))
    
    fluidRow(
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #1E3A8A; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", formatC(n_cells, format = "d", big.mark = ",")),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Inspected Cells (N)"))),
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #0D9488; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", formatC(n_feats, format = "d", big.mark = ",")),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Features (P)"))),
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #D97706; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", sprintf("%.1f%%", sparsity_val)),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Matrix Sparsity (S%)"))),
      column(3, div(class = "stat-card", style = "border: 1px solid #E2E8F0; border-left: 3px solid #6366F1; padding: 10px 12px; margin-bottom: 0; background: #FFFFFF; border-radius: 6px;",
                    div(class = "stat-number", style = "color: #0F172A; font-size: 1.35rem; font-family: monospace; font-weight: 700;", formatC(med_depth, format = "d", big.mark = ",")),
                    div(class = "stat-label", style = "font-size: 0.72rem; letter-spacing: 0.05em; color: #64748B;", "Median Library Depth (L~)")))
    )
  })

  output$table_matrix_inspect <- renderDT({
    res <- active_inspect_data()
    req(res, res$df)
    dt <- datatable(
      res$df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = TRUE,
        dom = "ftip"
      ),
      rownames = FALSE,
      class = "compact stripe hover border"
    )
    num_cols <- which(sapply(res$df, is.numeric))
    if (length(num_cols) > 0) {
      dt <- dt %>% formatStyle(num_cols, fontFamily = "monospace", textAlign = "right")
    }
    dt
  })

  output$download_inspect_data_csv <- downloadHandler(
    filename = function() {
      mode <- input$opt_inspect_view_mode %||% "view"
      paste0("scSimEval_inspect_", mode, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      res <- active_inspect_data()
      req(res, res$df)
      utils::write.csv(res$df, file, row.names = FALSE)
    }
  )

  output$ui_mset_preview_controls <- renderUI({
    req(rv$benchmark_df)
    if (!"Dataset" %in% colnames(rv$benchmark_df) || length(unique(rv$benchmark_df$Dataset)) <= 1) {
      return(NULL)
    }
    dsets <- unique(rv$benchmark_df$Dataset)
    div(
      style = "background: #F8FAFC; border: 1px solid #E2E8F0; border-radius: 6px; padding: 10px 14px; margin-bottom: 12px;",
      fluidRow(
        column(
          6,
          selectInput(
            "sel_preview_dataset", "View Scores For:",
            choices = c(
              "Cross-Dataset Average (Summary Across All Datasets)" = "avg",
              "Granular (All Datasets Combined)" = "all",
              setNames(dsets, paste0("Dataset: ", dsets))
            ),
            selected = "avg"
          )
        ),
        column(
          6,
          div(
            style = "margin-top: 25px; text-align: right;",
            downloadButton("download_preview_table_csv", "Export Table (CSV)", class = "btn btn-sm btn-outline-secondary", icon = icon("file-csv"))
          )
        )
      )
    )
  })

  output$download_preview_table_csv <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_preview_", Sys.Date(), ".csv") },
    content = function(file) {
      req(rv$benchmark_df)
      df <- rv$benchmark_df
      sel <- input$sel_preview_dataset %||% "avg"
      if (sel == "avg" && "Dataset" %in% colnames(df)) {
        has_norm <- "Score_Norm" %in% colnames(df)
        has_raw  <- "Score_Raw" %in% colnames(df)
        score_col <- if (has_raw) "Score_Raw" else if ("Score" %in% colnames(df)) "Score" else "Value"
        agg <- stats::aggregate(df[[score_col]], by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)
        colnames(agg)[4] <- "Average_Score"
        df <- agg
      } else if (sel != "all" && "Dataset" %in% colnames(df)) {
        df <- df[df$Dataset == sel, , drop = FALSE]
      }
      utils::write.csv(df, file, row.names = FALSE)
    }
  )

  output$table_active_data_preview <- renderDT({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    sel <- input$sel_preview_dataset %||% "avg"
    if ("Dataset" %in% colnames(df) && length(unique(df$Dataset)) > 1) {
      if (sel == "avg") {
        has_norm <- "Score_Norm" %in% colnames(df)
        has_raw  <- "Score_Raw" %in% colnames(df)
        score_col <- if (has_raw) "Score_Raw" else if ("Score" %in% colnames(df)) "Score" else "Value"
        agg <- stats::aggregate(df[[score_col]], by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)
        colnames(agg)[4] <- "Average_Score"
        df <- agg
      } else if (sel != "all") {
        df <- df[df$Dataset == sel, , drop = FALSE]
      }
    }
    datatable(
      head(df, 200),
      options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE),
      rownames = FALSE,
      class = "compact stripe hover border"
    ) %>% formatRound(columns = which(sapply(head(df, 200), is.numeric)), digits = 4)
  })

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # Tab 3: Comparative Bubble Matrix & Leaderboard (Minimalist & User-Friendly)
  # ----------------------------------------------------------------------------
  # Reactive list of available datasets for Comparative Bubble Matrix
  available_bubble_datasets <- reactive({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    if ("Dataset" %in% colnames(df)) {
      dsets <- unique(as.character(df$Dataset))
      dsets <- dsets[!is.na(dsets) & dsets != ""]
      if (length(dsets) > 0) return(dsets)
    }
    return("Dataset 1")
  })

  output$ui_bubble_dataset_selector <- renderUI({
    dsets <- available_bubble_datasets()
    choices <- list()
    if (length(dsets) > 1) {
      choices[[sprintf("All Datasets (Average Across %d Datasets)", length(dsets))]] <- "all_avg"
      for (d in dsets) {
        choices[[d]] <- d
      }
    } else if (length(dsets) == 1) {
      choices[["All Datasets (Average)"]] <- "all_avg"
      choices[[dsets[1]]] <- dsets[1]
    } else {
      choices[["All Datasets (Average)"]] <- "all_avg"
    }
    
    current_sel <- input$sel_bubble_dataset_view %||% "all_avg"
    if (!current_sel %in% unlist(choices)) current_sel <- "all_avg"
    
    div(
      class = "mb-2",
      tags$label(
        class = "control-label fw-bold mb-1",
        style = "font-size: 0.86rem; color: #1E293B;",
        icon("layer-group", class = "me-1 text-primary"),
        "Dataset Selection:"
      ),
      selectInput(
        "sel_bubble_dataset_view", NULL,
        choices = choices,
        selected = current_sel,
        width = "100%"
      )
    )
  })

  output$ui_bubble_modality_selector <- renderUI({
    req(rv$benchmark_df)
    has_multi <- "Modality" %in% colnames(rv$benchmark_df) &&
      (any(grepl("RNA", rv$benchmark_df$Modality, ignore.case = TRUE)) && any(grepl("ATAC", rv$benchmark_df$Modality, ignore.case = TRUE)))
    if (!has_multi) return(NULL)
    
    div(
      class = "mb-2",
      tags$label(
        class = "control-label fw-bold mb-1",
        style = "font-size: 0.86rem; color: #1E293B;",
        icon("dna", class = "me-1 text-primary"),
        "Assay Modality:"
      ),
      selectInput(
        "sel_bubble_modality", NULL,
        choices = c(
          "All Modalities (Combined Multiomics)" = "all",
          "scRNA-seq (Transcriptome)" = "RNA",
          "scATAC-seq (Chromatin Accessibility)" = "ATAC",
          "Cross-Modal Coupling (Joint Multiomics)" = "Joint"
        ),
        selected = input$sel_bubble_modality %||% (input$sel_viz_modality %||% "all"),
        width = "100%"
      )
    )
  })

  observeEvent(input$sel_bubble_modality, {
    if (!is.null(input$sel_viz_modality) && input$sel_viz_modality != input$sel_bubble_modality) {
      updateSelectInput(session, "sel_viz_modality", selected = input$sel_bubble_modality)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$sel_viz_modality, {
    if (!is.null(input$sel_bubble_modality) && input$sel_bubble_modality != input$sel_viz_modality) {
      updateSelectInput(session, "sel_bubble_modality", selected = input$sel_viz_modality)
    }
  }, ignoreInit = TRUE)

  output$ui_bubble_dataset_pills <- renderUI({
    dsets <- available_bubble_datasets()
    if (length(dsets) <= 1) return(NULL)
    
    current_sel <- input$sel_bubble_dataset_view %||% "all_avg"
    
    btn_all <- actionButton(
      "btn_view_dataset_all",
      label = tagList(icon("chart-pie", class = "me-1"), sprintf("All Datasets (Average Across %d Datasets)", length(dsets))),
      class = paste0("btn btn-sm me-1 mb-1 ", if (identical(current_sel, "all_avg")) "btn-primary active fw-bold" else "btn-outline-secondary"),
      onclick = "Shiny.setInputValue('sel_bubble_dataset_view', 'all_avg');"
    )
    
    dset_btns <- lapply(seq_along(dsets), function(i) {
      d <- dsets[i]
      actionButton(
        paste0("btn_view_dataset_", i),
        label = tagList(icon("database", class = "me-1"), d),
        class = paste0("btn btn-sm me-1 mb-1 ", if (identical(current_sel, d)) "btn-primary active fw-bold" else "btn-outline-secondary"),
        onclick = sprintf("Shiny.setInputValue('sel_bubble_dataset_view', '%s');", d)
      )
    })
    
    div(
      class = "p-2 mb-3 bg-light border rounded d-flex align-items-center flex-wrap gap-1",
      style = "border-left: 3px solid #2563EB !important;",
      div(class = "fw-bold text-dark me-2", style = "font-size: 0.84rem;", icon("table-columns", class = "me-1 text-primary"), "Dataset Tabs:"),
      btn_all,
      dset_btns
    )
  })

  observeEvent(input$btn_view_dataset_all, {
    updateSelectInput(session, "sel_bubble_dataset_view", selected = "all_avg")
  })

  output$ui_bubble_active_badge <- renderUI({
    df <- filtered_bubble_data()
    if (is.null(df) || nrow(df) == 0) {
      return(span(class = "badge bg-secondary", "No Metrics Selected"))
    }
    n_m <- length(unique(df$Metric))
    n_s <- length(unique(df$Method))
    d_sel <- input$sel_bubble_dataset_view %||% "all_avg"
    n_d <- if (!is.null(rv$benchmark_df) && "Dataset" %in% colnames(rv$benchmark_df)) length(unique(rv$benchmark_df$Dataset)) else 1
    d_tag <- if (identical(d_sel, "all_avg") && n_d > 1) sprintf("All %d Datasets (Average)", n_d) else if (identical(d_sel, "all_avg")) "All Datasets" else d_sel
    mode <- input$rad_bubble_filter_mode %||% "all"
    lbl <- if (identical(mode, "custom")) "Custom" else "All 62"
    mod_tag <- input$sel_bubble_modality %||% (input$sel_viz_modality %||% "all")
    mod_lbl <- if (!identical(mod_tag, "all")) sprintf(" | Modality: %s", switch(mod_tag, RNA = "scRNA-seq", ATAC = "scATAC-seq", Joint = "Cross-Modal", mod_tag)) else ""
    span(class = "badge bg-primary", style = "font-size: 0.8rem; font-weight: 500;",
         sprintf("%s%s | %s: %d Metrics x %d Simulators", d_tag, mod_lbl, lbl, n_m, n_s))
  })

  output$ui_bubble_method_picker <- renderUI({
    req(rv$methods)
    div(
      div(
        class = "d-flex justify-content-between align-items-center mb-1",
        tags$b("Simulators to Compare:"),
        div(
          actionLink("btn_bubble_methods_all", "All", style = "font-size: 0.78rem; margin-right: 6px; text-decoration: none; font-weight: 600;"),
          actionLink("btn_bubble_methods_none", "Clear", style = "font-size: 0.78rem; text-decoration: none; color: #dc3545;")
        )
      ),
      checkboxGroupInput(
        "sel_bubble_methods", NULL,
        choices = rv$methods,
        selected = rv$methods
      )
    )
  })

  observeEvent(input$btn_bubble_methods_all, {
    req(rv$methods)
    updateCheckboxGroupInput(session, "sel_bubble_methods", selected = rv$methods)
  })

  observeEvent(input$btn_bubble_methods_none, {
    updateCheckboxGroupInput(session, "sel_bubble_methods", selected = character(0))
  })

  output$ui_custom_category_picker <- renderUI({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    cats <- if ("Category" %in% colnames(df)) {
      avail <- unique(df$Category)
      avail <- avail[!is.na(avail) & trimws(avail) != ""]
      canon_order <- c(
        "(I) Distributional Properties",
        "(II) Correlations & Zero-Inflation",
        "(III) Cellular Structure & Concordance",
        "(IV) Batch Effects & Confounder Mixing",
        "(V) Biological Signal & Downstream Fidelity",
        "(VI) Trajectory & Lineage Dynamics",
        "(VII) Cross-Modal Coupling & Modularity",
        "(VIII) Computational Scalability"
      )
      intersect(canon_order, avail)
    } else {
      character(0)
    }
    choices <- c("All Categories (All Measures)" = "all", cats)
    selectInput("sel_custom_cat", "1. Filter by Category:", choices = choices, selected = "all")
  })

  output$ui_custom_metrics_checklist <- renderUI({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    cat_chosen <- input$sel_custom_cat %||% "all"
    
    sub_df <- if (identical(cat_chosen, "all") || !"Category" %in% colnames(df)) {
      df
    } else {
      df[df$Category == cat_chosen, , drop = FALSE]
    }
    
    avail_metrics <- unique(as.character(sub_df$Metric))
    avail_metrics <- avail_metrics[!is.na(avail_metrics) & trimws(avail_metrics) != ""]
    
    checkboxGroupInput(
      "chk_custom_metrics", NULL,
      choices = avail_metrics,
      selected = avail_metrics
    )
  })

  output$txt_custom_metrics_count <- renderText({
    n_sel <- length(input$chk_custom_metrics)
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    cat_chosen <- input$sel_custom_cat %||% "all"
    sub_df <- if (identical(cat_chosen, "all") || !"Category" %in% colnames(df)) df else df[df$Category == cat_chosen, , drop = FALSE]
    n_cat <- length(unique(sub_df$Metric))
    sprintf("%d / %d metrics selected", n_sel, n_cat)
  })

  observeEvent(input$btn_custom_metrics_all, {
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    cat_chosen <- input$sel_custom_cat %||% "all"
    sub_df <- if (identical(cat_chosen, "all") || !"Category" %in% colnames(df)) df else df[df$Category == cat_chosen, , drop = FALSE]
    avail_metrics <- unique(as.character(sub_df$Metric))
    updateCheckboxGroupInput(session, "chk_custom_metrics", selected = avail_metrics)
  })

  observeEvent(input$btn_custom_metrics_none, {
    updateCheckboxGroupInput(session, "chk_custom_metrics", selected = character(0))
  })

  filtered_bubble_data <- reactive({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    if (!"Dataset" %in% colnames(df)) {
      df$Dataset <- "Dataset 1"
    }
    
    # 1. Dataset Scope Filter & Cross-Dataset Averaging
    dataset_choice <- input$sel_bubble_dataset_view %||% "all_avg"
    all_dsets <- unique(as.character(df$Dataset))
    all_dsets <- all_dsets[!is.na(all_dsets) & all_dsets != ""]
    
    if (identical(dataset_choice, "all_avg") && length(all_dsets) > 1) {
      score_col <- if ("Score" %in% colnames(df)) "Score" else "Value"
      
      agg_score <- stats::aggregate(
        df[[score_col]],
        by = list(Method = df$Method, Category = df$Category, Metric = df$Metric),
        FUN = mean,
        na.rm = TRUE
      )
      colnames(agg_score)[4] <- score_col
      agg_score$Dataset <- sprintf("All Datasets (Average of %d Datasets)", length(all_dsets))
      
      if ("Score_Raw" %in% colnames(df)) {
        agg_raw <- stats::aggregate(df$Score_Raw, by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)
        agg_score$Score_Raw <- agg_raw[[4]]
      }
      if ("Score_Norm" %in% colnames(df)) {
        agg_norm <- stats::aggregate(df$Score_Norm, by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)
        agg_score$Score_Norm <- agg_norm[[4]]
      }
      df <- agg_score
    } else if (!identical(dataset_choice, "all_avg") && dataset_choice %in% all_dsets) {
      df <- df[df$Dataset == dataset_choice, , drop = FALSE]
    }
    
    # 2. Modality Filter
    mod_choice <- input$sel_bubble_modality %||% (input$sel_viz_modality %||% "all")
    if ("Modality" %in% colnames(df) && !identical(mod_choice, "all")) {
      if (identical(mod_choice, "RNA")) {
        df <- df[df$Modality %in% c("RNA", "Computation"), , drop = FALSE]
      } else if (identical(mod_choice, "ATAC")) {
        df <- df[df$Modality %in% c("ATAC", "Computation"), , drop = FALSE]
      } else if (identical(mod_choice, "Joint")) {
        df <- df[df$Modality %in% c("Joint (Cross-Modal)", "Computation"), , drop = FALSE]
      }
    }
    
    # 3. Simulator Methods Filter
    if (!is.null(input$sel_bubble_methods)) {
      if (length(input$sel_bubble_methods) == 0) return(NULL)
      df <- df[df$Method %in% input$sel_bubble_methods, , drop = FALSE]
    }
    
    # 3. Mode Filter: 'all' vs 'custom'
    mode <- input$rad_bubble_filter_mode %||% "all"
    
    if (identical(mode, "custom")) {
      cat_chosen <- input$sel_custom_cat %||% "all"
      if (!identical(cat_chosen, "all") && "Category" %in% colnames(df)) {
        df <- df[df$Category == cat_chosen, , drop = FALSE]
      }
      
      if (!is.null(input$chk_custom_metrics)) {
        if (length(input$chk_custom_metrics) == 0) return(NULL)
        if ("Metric" %in% colnames(df)) {
          df <- df[df$Metric %in% input$chk_custom_metrics, , drop = FALSE]
        }
      }
    }
    
    if (nrow(df) == 0) return(NULL)
    if ("Metric" %in% colnames(df)) df <- df[!is.na(df$Metric) & df$Metric != "" & df$Metric != "NA", , drop = FALSE]
    if ("Category" %in% colnames(df)) df <- df[!is.na(df$Category) & df$Category != "" & df$Category != "NA", , drop = FALSE]
    score_c <- if ("Score" %in% colnames(df)) "Score" else "Value"
    if (score_c %in% colnames(df)) df <- df[!is.na(df[[score_c]]), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df
  })
  bubble_plot_reactive <- reactive({
    req(rv_plot_triggers$bubble > 0)
    df <- filtered_bubble_data()
    if (is.null(df) || nrow(df) == 0) {
      p_empty <- ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 1, y = 1,
          label = "No metrics or simulators selected.\nPlease select at least one simulator or check metrics in Custom Selection.",
          size = 5.2, color = "#64748B", fontface = "italic"
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.background = ggplot2::element_rect(fill = "#F8FAFC", color = "#E2E8F0"))
      return(p_empty)
    }
    
    withProgress(message = "Rendering Comparative Bubble Matrix", detail = "Normalizing 62 metrics...", {
      df <- standardize_benchmark_categories(df)
      incProgress(0.5, detail = "Generating ggplot bubble matrix...")
      p <- plot_benchmark_bubble_matrix(
        data              = df,
        base_size         = 11,
        compact_strips    = TRUE,
        show_missing_dots = FALSE,
        normalize_scores  = TRUE
      )
      
      n_metrics <- length(unique(df$Metric))
      n_cats <- length(unique(df$Category))
      n_sims <- length(unique(df$Method))
      mode <- input$rad_bubble_filter_mode %||% "all"
      dataset_choice <- input$sel_bubble_dataset_view %||% "all_avg"
      all_dsets <- if (!is.null(rv$benchmark_df) && "Dataset" %in% colnames(rv$benchmark_df)) unique(as.character(rv$benchmark_df$Dataset)) else "Dataset 1"
      is_avg <- identical(dataset_choice, "all_avg")
      n_d <- length(all_dsets)
      
      # Clean publication-ready title & subtitle
      show_title <- if (!is.null(input$chk_bubble_show_title)) isTRUE(input$chk_bubble_show_title) else TRUE
      custom_title <- input$txt_bubble_custom_title
      
      main_title <- if (!show_title) {
        NULL
      } else if (!is.null(custom_title) && trimws(custom_title) != "") {
        trimws(custom_title)
      } else if (is_avg && n_d > 1) {
        "Benchmark Evaluation of Single-Cell Simulators"
      } else if (!is_avg) {
        sprintf("Comparative Evaluation of Single-Cell Simulators (%s)", dataset_choice)
      } else {
        "Comparative Evaluation of Single-Cell Simulators"
      }
      
      sub_title <- if (!show_title) {
        NULL
      } else if (is_avg && n_d > 1) {
        sprintf("Standardized fidelity scores [0, 1] averaged across %d benchmark datasets", n_d)
      } else if (!is_avg) {
        sprintf("Standardized fidelity scores [0, 1] on %s", dataset_choice)
      } else {
        "Standardized fidelity scores [0, 1] across evaluation tiers"
      }
      
      caption_text <- if (!show_title) NULL else paste0(
        "Circle: standard performance  |  Bold circle outline: top-performing simulator for the metric  |  Higher score = closer concordance with empirical reference."
      )
      
      p + ggplot2::labs(
        title = main_title,
        subtitle = sub_title,
        caption = caption_text
      )
    })
  })

  output$ui_bubble_plot_render <- renderUI({
    if (rv_plot_triggers$bubble == 0) {
      return(render_plot_placeholder(
        title = "Comparative Bubble Matrix Ready to Plot",
        desc = "Choose your simulators and metric categories on the left, then click below to confirm and render the flagship comparative bubble matrix.",
        button_id = "btn_plot_bubble_inline",
        btn_label = "Confirm & Generate Bubble Matrix",
        icon_name = "circle-dot"
      ))
    }
    df <- filtered_bubble_data()
    n_m <- if (!is.null(df)) length(unique(df$Metric)) else 20
    default_h <- max(380, min(1200, 160 + n_m * 22))
    default_w <- if (n_m <= 10) 1400 else 2200
    h_val <- if (!is.null(input$sld_bubble_height)) input$sld_bubble_height else default_h
    w_val <- if (!is.null(input$sld_bubble_width)) input$sld_bubble_width else default_w
    plotOutput("plot_bubble_matrix", width = paste0(w_val, "px"), height = paste0(h_val, "px"))
  })

  output$plot_bubble_matrix <- renderPlot({
    bubble_plot_reactive()
  })

  leaderboard_reactive <- reactive({
    df <- filtered_bubble_data()
    if (is.null(df) || nrow(df) == 0) {
      return(data.frame(
        Overall_Rank = integer(0),
        Method = character(0),
        Average_Fidelity = character(0),
        Fidelity_Score = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    compute_method_leaderboard(df)
  })

  output$table_leaderboard_dt <- renderDT({
    lb <- leaderboard_reactive()
    if (nrow(lb) == 0) {
      return(datatable(data.frame(Notice = "No data selected for ranking. Please select simulators or check metrics from the sidebar.")))
    }
    display_df <- lb
    colnames(display_df) <- c("Overall Rank", "Method", "Average Fidelity", "Fidelity Score")
    
    datatable(
      display_df,
      options = list(
        dom = "t",
        pageLength = 25,
        ordering = FALSE,
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      ),
      rownames = FALSE,
      class = "compact stripe hover border"
    ) %>%
      formatStyle(
        "Overall Rank",
        fontWeight = "bold",
        backgroundColor = styleEqual(
          c(1, 2, 3),
          c("#FFF9DB", "#F1F3F5", "#FFF4E6")
        )
      ) %>%
      formatStyle(
        "Average Fidelity",
        fontWeight = "bold",
        color = "#1D72B8"
      ) %>%
      formatStyle(
        "Fidelity Score",
        fontFamily = "monospace",
        fontWeight = "bold"
      )
  })

  output$download_leaderboard_csv <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_leaderboard_", Sys.Date(), ".csv") },
    content = function(file) {
      lb <- leaderboard_reactive()
      colnames(lb) <- c("Overall Rank", "Method", "Average Fidelity", "Fidelity Score")
      utils::write.csv(lb, file, row.names = FALSE)
    }
  )

  output$download_leaderboard_csv_side <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_leaderboard_", Sys.Date(), ".csv") },
    content = function(file) {
      lb <- leaderboard_reactive()
      colnames(lb) <- c("Overall Rank", "Method", "Average Fidelity", "Fidelity Score")
      utils::write.csv(lb, file, row.names = FALSE)
    }
  )

  output$download_bubble_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_bubble_matrix_", make.names(input$sel_bubble_dataset_view %||% "All_Datasets_Average"), "_", Sys.Date(), ".jpeg") },
    content = function(file) {
      df <- filtered_bubble_data()
      n_m <- if (!is.null(df)) length(unique(df$Metric)) else 20
      h_in <- max(4.5, min(14, 2.5 + n_m * 0.16))
      w_in <- if (!is.null(input$sld_bubble_width)) max(10, input$sld_bubble_width / 95) else 23
      export_single_jpeg(file, bubble_plot_reactive(), width = w_in, height = h_in, dpi = 600)
    }
  )

  output$download_bubble_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_bubble_matrix_", make.names(input$sel_bubble_dataset_view %||% "All_Datasets_Average"), "_", Sys.Date(), ".pdf") },
    content = function(file) {
      df <- filtered_bubble_data()
      n_m <- if (!is.null(df)) length(unique(df$Metric)) else 20
      h_in <- max(4.5, min(14, 2.5 + n_m * 0.16))
      w_in <- if (!is.null(input$sld_bubble_width)) max(10, input$sld_bubble_width / 95) else 23
      grDevices::pdf(file, width = w_in, height = h_in)
      print(bubble_plot_reactive())
      grDevices::dev.off()
    }
  )

  # Helper: Generate Comprehensive File Path Configuration Block for Downloaded R Scripts
  get_workflow_file_paths_block <- function(target_analysis = "benchmark") {
    lines <- c(
      "",
      "# ==============================================================================",
      "# ==============================================================================",
      "# FILE PATH CONFIGURATION & HPC CLUSTER RUN GUIDE",
      "# ==============================================================================",
      "# Instructions for Non-Coding Users & HPC / Cluster Deployment:",
      "#   1. Copy this R script and all your data files into your working folder on",
      "#      your HPC cluster, cloud instance, or server.",
      "#   2. Update the file path variables below to point to the exact locations on your system.",
      "#   3. Run the script directly from your terminal:",
      "#        Rscript <this_script_name>.R",
      "#      Or submit as a SLURM cluster job:",
      "#        sbatch --cpus-per-task=8 --mem=32G --time=02:00:00 --wrap=\"Rscript <this_script_name>.R\"",
      "# ==============================================================================",
      "",
      "# 1. Primary Benchmark Results Table (Exported from scSimEval):",
      "# Point this variable to your exported CSV or RDS results file from the Shiny app:",
      "BENCHMARK_RESULTS_FILE <- 'scSimEval_benchmark_results.csv' # or 'scSimEval_benchmark_results.rds'",
      ""
    )
    
    is_multi <- !is.null(rv$benchmark_df) && "Modality" %in% colnames(rv$benchmark_df) &&
      (any(grepl("RNA", rv$benchmark_df$Modality, ignore.case = TRUE)) && any(grepl("ATAC", rv$benchmark_df$Modality, ignore.case = TRUE)))
    
    if (is_multi) {
      lines <- c(lines,
        "# 2. Input Single-Cell Multiomics Count Matrices & Annotations:",
        "# Update paths to your empirical reference and simulated datasets on HPC/cloud:",
        "REF_RNA_FILE    <- 'reference_rna.csv'   # (Supported: .csv, .tsv, .rds, .mtx, .xlsx)",
        "REF_ATAC_FILE   <- 'reference_atac.csv'  # (Supported: .csv, .tsv, .rds, .mtx, .xlsx)",
        "REF_DATA_FILE   <- REF_RNA_FILE          # Alias for Distribution QC workflow",
        "REF_MATRIX_FILE <- REF_RNA_FILE          # Alias for Cell Embeddings workflow",
        "CELL_TYPES_FILE <- 'cell_types.csv'      # (Optional: cell barcode to cell type mapping)",
        "BATCH_FILE      <- 'batches.csv'         # (Optional: cell barcode to batch mapping)",
        "",
        "# Simulated Datasets for Evaluated Simulators:"
      )
      sim_names <- if (!is.null(rv$methods) && length(rv$methods) > 0) rv$methods else c("scMultiSim", "scDesign3", "MOSim", "scCross")
      lines <- c(lines, "SIMULATOR_FILES <- list(")
      for (i in seq_along(sim_names)) {
        m <- sim_names[i]
        comma <- if (i < length(sim_names)) "," else ""
        lines <- c(lines, sprintf("  '%s' = list(rna = 'sim_%s_rna.csv', atac = 'sim_%s_atac.csv')%s", m, m, m, comma))
      }
      lines <- c(lines, ")",
        "SIM_DATA_FILE   <- 'sim_scMultiSim_rna.csv' # Primary simulator for Distribution QC",
        ""
      )
    } else {
      lines <- c(lines,
        "# 2. Input Single-Cell Count Matrix & Annotations:",
        "# Update paths to your empirical reference and simulated datasets on HPC/cloud:",
        "REF_MATRIX_FILE <- 'reference_counts.csv' # (Supported: .csv, .tsv, .rds, .mtx, .xlsx)",
        "REF_DATA_FILE   <- REF_MATRIX_FILE        # Alias for Distribution QC workflow",
        "CELL_TYPES_FILE <- 'cell_types.csv'       # (Optional: cell barcode to cell type mapping)",
        "BATCH_FILE      <- 'batches.csv'          # (Optional: cell barcode to batch mapping)",
        "",
        "# Simulated Count Matrices for Evaluated Simulators:"
      )
      sim_names <- if (!is.null(rv$methods) && length(rv$methods) > 0) rv$methods else c("scDesign3", "Splatter", "SCRIP", "SymSim")
      lines <- c(lines, "SIMULATOR_FILES <- list(")
      for (i in seq_along(sim_names)) {
        m <- sim_names[i]
        comma <- if (i < length(sim_names)) "," else ""
        lines <- c(lines, sprintf("  '%s' = 'sim_%s_counts.csv'%s", m, m, comma))
      }
      lines <- c(lines, ")",
        "SIM_DATA_FILE   <- 'sim_scDesign3_counts.csv' # Primary simulator for Distribution QC",
        ""
      )
    }
    
    lines <- c(lines,
      "# 3. HPC Parallel Acceleration & Output Configuration:",
      "NUM_THREADS <- 4L                   # Number of CPU cores for parallel metric evaluation",
      "OUTPUT_DIR  <- './scSimEval_output' # Directory where figures and results tables will be saved",
      "if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)",
      "# ==============================================================================",
      ""
    )
    return(paste(lines, collapse = "\n"))
  }

  # Generator for Standalone Reproducible Analytical Workflow Script for Comparative Bubble Matrix
  generate_bubble_plot_script <- function() {
    date_str <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    df_raw <- rv$benchmark_df
    if (is.null(df_raw) || nrow(df_raw) == 0) {
      return("#!/usr/bin/env Rscript\n# scSimEval Comparative Bubble Matrix Analysis Workflow\n# Note: No benchmark data currently loaded.\n")
    }
    
    dataset_choice <- input$sel_bubble_dataset_view %||% "all_avg"
    all_dsets <- if ("Dataset" %in% colnames(df_raw)) unique(as.character(df_raw$Dataset)) else "Dataset 1"
    all_dsets <- all_dsets[!is.na(all_dsets) & all_dsets != ""]
    
    selected_methods <- input$sel_bubble_methods
    if (is.null(selected_methods) || length(selected_methods) == 0) {
      selected_methods <- unique(as.character(df_raw$Method))
    }
    
    filter_mode <- input$rad_bubble_filter_mode %||% "all"
    selected_cat <- input$sel_custom_cat %||% "all"
    selected_metrics <- input$chk_custom_metrics
    
    show_title <- if (!is.null(input$chk_bubble_show_title)) isTRUE(input$chk_bubble_show_title) else TRUE
    custom_title <- input$txt_bubble_custom_title %||% ""
    
    w_px <- input$sld_bubble_width %||% 2200
    h_px <- input$sld_bubble_height %||% 680
    w_in <- round(max(10, w_px / 95), 2)
    h_in <- round(max(4.5, min(14, 2.5 + 35 * 0.16)), 2)
    
    methods_code <- paste0("c(", paste(sprintf("\"%s\"", selected_methods), collapse = ", "), ")")
    metrics_code <- if (!is.null(selected_metrics) && length(selected_metrics) > 0) {
      paste0("c(", paste(sprintf("\"%s\"", selected_metrics), collapse = ", "), ")")
    } else {
      "NULL"
    }
    
    out <- c(
      "#!/usr/bin/env Rscript",
      "# ==============================================================================",
      "# scSimEval: Comparative Benchmark Bubble Matrix Analysis Workflow",
      sprintf("# Generated on: %s", date_str),
      sprintf("# Dataset Scope: %s", dataset_choice),
      sprintf("# Simulators: %s", paste(selected_methods, collapse = ", ")),
      sprintf("# Filter Mode: %s", filter_mode),
      "# ==============================================================================",
      "# Analytical Workflow for evaluating single-cell simulation methods across",
      "# multi-tiered benchmark measures with standardized directional scoring.",
      "# ==============================================================================",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW STEP 1: Environment Setup & Package Dependencies",
      "# ------------------------------------------------------------------------------",
      "suppressPackageStartupMessages({",
      "  library(scSimEval)",
      "  library(ggplot2)",
      "})",
      "",
      "cat('[scSimEval Workflow] Step 1: Environment initialized successfully.\n')",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW STEP 2: Ingest Benchmark Evaluation Data",
      "# ------------------------------------------------------------------------------",
      "# Check if benchmark results file exists (see FILE PATH CONFIGURATION at end of script)",
      "if (!exists('BENCHMARK_RESULTS_FILE')) BENCHMARK_RESULTS_FILE <- 'scSimEval_benchmark_results.csv'",
      "",
      "if (file.exists(BENCHMARK_RESULTS_FILE)) {",
      "  cat(sprintf('[scSimEval Workflow] Step 2: Loading benchmark evaluation data from: %s\n', BENCHMARK_RESULTS_FILE))",
      "  if (endsWith(tolower(BENCHMARK_RESULTS_FILE), '.rds')) {",
      "    obj <- readRDS(BENCHMARK_RESULTS_FILE)",
      "    benchmark_data <- if (is.data.frame(obj)) obj else (obj$benchmark_summary_table %||% obj$benchmark_df)",
      "  } else {",
      "    benchmark_data <- utils::read.csv(BENCHMARK_RESULTS_FILE, stringsAsFactors = FALSE, check.names = FALSE)",
      "  }",
      "} else {",
      "  cat(sprintf('[scSimEval Workflow] Step 2: BENCHMARK_RESULTS_FILE (\"%s\") not found on disk.\n', BENCHMARK_RESULTS_FILE))",
      "  cat('  Please set the correct path in the FILE PATH CONFIGURATION section at the bottom of this script.\n')",
      "  cat('  Loading package demo benchmark database for execution demonstration...\n')",
      "  pth <- system.file('shiny/scSimEvalApp/data/demo_benchmark_data.rds', package = 'scSimEval')",
      "  if (file.exists(pth)) {",
      "    demo_obj <- readRDS(pth)",
      "    benchmark_data <- demo_obj$benchmark_summary_table %||% demo_obj$benchmark_df",
      "  } else {",
      "    stop(sprintf('File \"%s\" not found. Please set the correct path in FILE PATH CONFIGURATION at the end of the script.', BENCHMARK_RESULTS_FILE))",
      "  }",
      "}",
      "",
      "cat(sprintf('[scSimEval Workflow] Ingested %d benchmark observations across %d simulators.\n', nrow(benchmark_data), length(unique(benchmark_data$Method)))),",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW STEP 3: Multi-Dataset Harmonization & Cross-Dataset Aggregation",
      "# ------------------------------------------------------------------------------",
      sprintf("dataset_scope <- '%s'", dataset_choice),
      "df <- benchmark_data",
      "if (!'Dataset' %in% colnames(df)) df$Dataset <- 'Dataset 1'",
      "all_dsets <- unique(as.character(df$Dataset))",
      "all_dsets <- all_dsets[!is.na(all_dsets) & all_dsets != '']",
      "",
      "if (identical(dataset_scope, 'all_avg') && length(all_dsets) > 1) {",
      "  score_col <- if ('Score' %in% colnames(df)) 'Score' else 'Value'",
      "  cat(sprintf('[scSimEval Workflow] Step 3: Computing cross-dataset arithmetic average across %d datasets...\n', length(all_dsets)))",
      "  agg_score <- stats::aggregate(df[[score_col]], by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)",
      "  colnames(agg_score)[4] <- score_col",
      "  agg_score$Dataset <- sprintf('All Datasets (Average of %d Datasets)', length(all_dsets))",
      "  if ('Score_Raw' %in% colnames(df)) agg_score$Score_Raw <- stats::aggregate(df$Score_Raw, by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)[[4]]",
      "  if ('Score_Norm' %in% colnames(df)) agg_score$Score_Norm <- stats::aggregate(df$Score_Norm, by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)[[4]]",
      "  df <- agg_score",
      "} else if (!identical(dataset_scope, 'all_avg') && dataset_scope %in% all_dsets) {",
      "  cat(sprintf('[scSimEval Workflow] Step 3: Subsetting strictly to dataset: %s\n', dataset_scope))",
      "  df <- df[df$Dataset == dataset_scope, , drop = FALSE]",
      "}",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW STEP 4: Analytical Subsetting (Simulators & Evaluation Categories)",
      "# ------------------------------------------------------------------------------",
      sprintf("selected_sims <- %s", methods_code),
      sprintf("filter_mode   <- '%s'", filter_mode),
      sprintf("selected_cat  <- '%s'", selected_cat),
      sprintf("selected_mets <- %s", metrics_code),
      "",
      "if (!is.null(selected_sims) && length(selected_sims) > 0) {",
      "  df <- df[df$Method %in% selected_sims, , drop = FALSE]",
      "}",
      "if (identical(filter_mode, 'custom')) {",
      "  if (!identical(selected_cat, 'all') && 'Category' %in% colnames(df)) {",
      "    df <- df[df$Category == selected_cat, , drop = FALSE]",
      "  }",
      "  if (!is.null(selected_mets) && length(selected_mets) > 0 && 'Metric' %in% colnames(df)) {",
      "    df <- df[df$Metric %in% selected_mets, , drop = FALSE]",
      "  }",
      "}",
      "",
      "# Standardize canonical category names",
      "if (exists('standardize_benchmark_categories', where = asNamespace('scSimEval'), mode = 'function')) {",
      "  df <- scSimEval::standardize_benchmark_categories(df)",
      "}",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW STEP 5: High-Density Comparative Bubble Matrix Construction",
      "# ------------------------------------------------------------------------------",
      "cat('[scSimEval Workflow] Step 5: Rendering publication comparative bubble matrix...\n')",
      "p <- plot_benchmark_bubble_matrix(",
      "  data              = df,",
      "  base_size         = 11,",
      "  compact_strips    = TRUE,",
      "  show_missing_dots = FALSE,",
      "  normalize_scores  = TRUE",
      ")",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW STEP 6: Publication Labels & Annotation",
      "# ------------------------------------------------------------------------------",
      sprintf("show_title   <- %s", if (show_title) "TRUE" else "FALSE"),
      sprintf("custom_title <- '%s'", custom_title),
      "if (!show_title) {",
      "  p <- p + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)",
      "} else {",
      "  main_title <- if (nzchar(custom_title)) custom_title else if (identical(dataset_scope, 'all_avg') && length(all_dsets) > 1) 'Benchmark Evaluation of Single-Cell Simulators' else 'Comparative Evaluation of Single-Cell Simulators'",
      "  sub_title  <- if (identical(dataset_scope, 'all_avg') && length(all_dsets) > 1) sprintf('Standardized fidelity scores [0, 1] averaged across %d benchmark datasets', length(all_dsets)) else sprintf('Standardized fidelity scores [0, 1] across evaluation tiers')",
      "  caption_text <- 'Circle: standard | Bold circle outline: top-performing simulator for the metric | Higher score = closer concordance with empirical reference.'",
      "  p <- p + ggplot2::labs(title = main_title, subtitle = sub_title, caption = caption_text)",
      "}",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW STEP 7: Multi-Format High-Resolution Export (Vector PDF & 600 DPI JPEG)",
      "# ------------------------------------------------------------------------------",
      sprintf("export_width  <- %.2f", w_in),
      sprintf("export_height <- %.2f", h_in),
      "pdf_file  <- 'scSimEval_comparative_bubble_matrix.pdf'",
      "jpeg_file <- 'scSimEval_comparative_bubble_matrix.jpeg'",
      "",
      "ggplot2::ggsave(pdf_file, plot = p, width = export_width, height = export_height, device = 'pdf')",
      "ggplot2::ggsave(jpeg_file, plot = p, width = export_width, height = export_height, dpi = 600)",
      "",
      "cat(sprintf('[scSimEval Workflow] Complete! Figures successfully saved:\n  - Vector PDF: %s\n  - Publication JPEG (600 DPI): %s\n', pdf_file, jpeg_file))",
      "",
      get_workflow_file_paths_block("bubble")
    )
    return(paste(out, collapse = "\n"))
  }

  output$download_bubble_r_script <- downloadHandler(
    filename = function() {
      ds_tag <- if (identical(input$sel_bubble_dataset_view, "all_avg")) "all_datasets" else gsub("[^A-Za-z0-9_]", "_", input$sel_bubble_dataset_view %||% "dataset1")
      paste0("scSimEval_bubble_matrix_", ds_tag, "_", format(Sys.Date(), "%Y%m%d"), ".R")
    },
    content = function(file) {
      writeLines(generate_bubble_plot_script(), file)
    }
  )

  output$download_bubble_r_script_header <- downloadHandler(
    filename = function() {
      ds_tag <- if (identical(input$sel_bubble_dataset_view, "all_avg")) "all_datasets" else gsub("[^A-Za-z0-9_]", "_", input$sel_bubble_dataset_view %||% "dataset1")
      paste0("scSimEval_bubble_matrix_", ds_tag, "_", format(Sys.Date(), "%Y%m%d"), ".R")
    },
    content = function(file) {
      writeLines(generate_bubble_plot_script(), file)
    }
  )

  # ============================================================================
  # Tab 4: Diagnostic Visualizations (Multi-Dataset Support & R Script Export)
  # ============================================================================
  
  # Centralized multi-dataset filtering and cross-dataset arithmetic average helper
  get_viz_benchmark_data <- function(dataset_choice = NULL, modality_choice = NULL) {
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    if (!"Dataset" %in% colnames(df)) df$Dataset <- "Dataset 1"
    
    active_subtab <- input$viz_subtabs %||% "1. Evaluation Summary"
    tab_dset <- if (grepl("1\\.", active_subtab)) input$sel_sum_dataset else
                if (grepl("2\\.", active_subtab)) input$sel_dist_dataset else
                if (grepl("3\\.", active_subtab)) input$sel_scale_dataset else
                if (grepl("4\\.", active_subtab)) (if (identical(input$opt_box_view_mode, "category")) input$sel_box_dataset_cat else input$sel_box_dataset) else
                if (grepl("5\\.", active_subtab)) input$sel_heat_dataset else
                if (grepl("6\\.", active_subtab)) input$sel_pca_dataset else
                if (grepl("7\\.", active_subtab)) input$sel_mds_dataset else
                if (grepl("8\\.", active_subtab)) input$sel_emb_dataset else NULL
    choice <- dataset_choice %||% (tab_dset %||% (input$sel_viz_dataset_view %||% (input$sel_bubble_dataset_view %||% "all_avg")))
    all_dsets <- unique(as.character(df$Dataset))
    all_dsets <- all_dsets[!is.na(all_dsets) & all_dsets != ""]
    
    # 1. Modality Filter
    tab_mod <- if (grepl("1\\.", active_subtab)) input$sel_sum_modality else
               if (grepl("2\\.", active_subtab)) input$sel_dist_modality else
               if (grepl("3\\.", active_subtab)) input$sel_scale_modality else
               if (grepl("4\\.", active_subtab)) (if (identical(input$opt_box_view_mode, "category")) input$sel_box_modality_cat else input$sel_box_modality) else
               if (grepl("5\\.", active_subtab)) input$sel_heat_modality else
               if (grepl("6\\.", active_subtab)) input$sel_pca_modality else
               if (grepl("7\\.", active_subtab)) input$sel_mds_modality else
               if (grepl("8\\.", active_subtab)) input$sel_emb_modality else NULL
    mod_sel <- modality_choice %||% (tab_mod %||% (input$sel_viz_modality %||% (input$sel_bubble_modality %||% "all")))
    if ("Modality" %in% colnames(df) && !identical(mod_sel, "all")) {
      if (identical(mod_sel, "RNA") || identical(mod_sel, "rna")) {
        df <- df[grepl("RNA", df$Modality, ignore.case = TRUE) | df$Modality %in% c("Computation", "Total"), , drop = FALSE]
      } else if (identical(mod_sel, "ATAC") || identical(mod_sel, "atac")) {
        df <- df[grepl("ATAC", df$Modality, ignore.case = TRUE) | df$Modality %in% c("Computation", "Total"), , drop = FALSE]
      } else if (identical(mod_sel, "Joint")) {
        df <- df[grepl("Joint", df$Modality, ignore.case = TRUE) | df$Modality %in% c("Computation", "Total"), , drop = FALSE]
      }
    }
    
    if (identical(choice, "all_avg") && length(all_dsets) > 1) {
      score_col <- if ("Score" %in% colnames(df)) "Score" else "Value"
      agg_score <- stats::aggregate(
        df[[score_col]],
        by = list(Method = df$Method, Category = df$Category, Metric = df$Metric),
        FUN = mean,
        na.rm = TRUE
      )
      colnames(agg_score)[4] <- score_col
      agg_score$Dataset <- sprintf("All Datasets (Average of %d Datasets)", length(all_dsets))
      
      if ("Score_Raw" %in% colnames(df)) {
        agg_raw <- stats::aggregate(df$Score_Raw, by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)
        agg_score$Score_Raw <- agg_raw[[4]]
      }
      if ("Score_Norm" %in% colnames(df)) {
        agg_norm <- stats::aggregate(df$Score_Norm, by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)
        agg_score$Score_Norm <- agg_norm[[4]]
      }
      if ("Modality" %in% colnames(df)) {
        agg_mod <- stats::aggregate(df$Modality, by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = function(x) x[1])
        agg_score$Modality <- agg_mod[[4]]
      }
      return(agg_score)
    } else if (!identical(choice, "all_avg") && choice %in% all_dsets) {
      return(df[df$Dataset == choice, , drop = FALSE])
    }
    return(df)
  }
  
  # Global Visualizations Dataset Navigation Bar
  output$ui_viz_global_dataset_bar <- renderUI({
    dsets <- available_bubble_datasets()
    has_multi <- !is.null(rv$benchmark_df) && "Modality" %in% colnames(rv$benchmark_df) &&
      (any(grepl("RNA", rv$benchmark_df$Modality, ignore.case = TRUE)) && any(grepl("ATAC", rv$benchmark_df$Modality, ignore.case = TRUE)))
    if (length(dsets) <= 1 && !has_multi) return(NULL)
    
    current_sel <- input$sel_viz_dataset_view %||% (input$sel_bubble_dataset_view %||% "all_avg")
    
    btn_all <- if (length(dsets) > 1) {
      actionButton(
        "btn_viz_view_dataset_all",
        label = tagList(icon("chart-pie", class = "me-1"), sprintf("All Datasets (Average Across %d Datasets)", length(dsets))),
        class = paste0("btn btn-sm me-1 mb-1 ", if (identical(current_sel, "all_avg")) "btn-primary active fw-bold" else "btn-outline-secondary"),
        onclick = "Shiny.setInputValue('sel_viz_dataset_view', 'all_avg'); Shiny.setInputValue('sel_bubble_dataset_view', 'all_avg');"
      )
    } else NULL
    
    dset_btns <- if (length(dsets) > 1) {
      lapply(seq_along(dsets), function(i) {
        d <- dsets[i]
        actionButton(
          paste0("btn_viz_view_dataset_", i),
          label = tagList(icon("database", class = "me-1"), d),
          class = paste0("btn btn-sm me-1 mb-1 ", if (identical(current_sel, d)) "btn-primary active fw-bold" else "btn-outline-secondary"),
          onclick = sprintf("Shiny.setInputValue('sel_viz_dataset_view', '%s'); Shiny.setInputValue('sel_bubble_dataset_view', '%s');", d, d)
        )
      })
    } else NULL
    
    choices <- list()
    if (length(dsets) > 1) {
      choices[[sprintf("All Datasets (Average Across %d Datasets)", length(dsets))]] <- "all_avg"
      for (d in dsets) choices[[d]] <- d
    }
    
    mod_selector <- if (has_multi) {
      div(
        class = "d-flex align-items-center gap-1",
        span(class = "fw-bold text-dark", style = "font-size: 0.86rem;", icon("dna", class = "me-1 text-primary"), "Assay Modality:"),
        selectInput(
          "sel_viz_modality", NULL,
          choices = c(
            "All Modalities (Combined Multiomics)" = "all",
            "scRNA-seq (Transcriptome)" = "RNA",
            "scATAC-seq (Chromatin Accessibility)" = "ATAC",
            "Cross-Modal Coupling (Joint Multiomics)" = "Joint"
          ),
          selected = input$sel_viz_modality %||% "all",
          width = "230px"
        )
      )
    } else NULL
    
    div(
      class = "p-2 mb-3 bg-light border rounded d-flex justify-content-between align-items-center flex-wrap gap-2",
      style = "border-left: 4px solid #2563EB !important;",
      div(
        class = "d-flex align-items-center flex-wrap gap-1",
        if (length(dsets) > 1) div(class = "fw-bold text-dark me-2", style = "font-size: 0.86rem;", icon("layer-group", class = "me-1 text-primary"), "Active Dataset Scope:") else NULL,
        btn_all,
        dset_btns
      ),
      div(
        class = "d-flex align-items-center flex-wrap gap-2",
        mod_selector,
        if (length(dsets) > 1) span(class = "badge bg-primary px-2 py-1", style = "font-size: 0.80rem;",
             if (identical(current_sel, "all_avg")) sprintf("Viewing: Average (%d Datasets)", length(dsets)) else sprintf("Viewing: %s", current_sel)) else NULL,
        if (length(dsets) > 1) selectInput("sel_viz_dataset_view", NULL, choices = choices, selected = current_sel, width = "220px") else NULL
      )
    )
  })
  
  observeEvent(input$btn_viz_view_dataset_all, {
    updateSelectInput(session, "sel_viz_dataset_view", selected = "all_avg")
  })
  
  # Active Badge Renderers for all 8 Sub-panels
  lapply(c("sum", "dist", "scale", "box", "heat", "pca", "mds", "emb"), function(tag) {
    output[[paste0("ui_viz_active_badge_", tag)]] <- renderUI({
      d_sel <- switch(tag,
        sum   = input$sel_sum_dataset,
        dist  = input$sel_dist_dataset,
        scale = input$sel_scale_dataset,
        box   = input$sel_box_dataset,
        heat  = input$sel_heat_dataset,
        pca   = input$sel_pca_dataset,
        mds   = input$sel_mds_dataset,
        emb   = input$sel_emb_dataset,
        input$sel_viz_dataset_view %||% "all_avg"
      ) %||% (input$sel_viz_dataset_view %||% "all_avg")
      n_d <- if (!is.null(rv$benchmark_df) && "Dataset" %in% colnames(rv$benchmark_df)) length(unique(rv$benchmark_df$Dataset)) else 1
      label <- if (identical(d_sel, "all_avg") && n_d > 1) sprintf("All %d Datasets (Average)", n_d) else if (identical(d_sel, "all_avg")) "All Datasets" else d_sel
      
      mod_tag <- input$sel_viz_modality %||% (input$sel_bubble_modality %||% "all")
      mod_badge <- if (!identical(mod_tag, "all")) {
        span(class = "badge bg-info-subtle text-info border border-info-subtle ms-1", style = "font-size: 0.80rem;",
             icon("dna", class = "me-1"), switch(mod_tag, RNA = "scRNA-seq", ATAC = "scATAC-seq", Joint = "Cross-Modal", mod_tag))
      } else NULL
      
      tagList(
        span(class = "badge bg-primary-subtle text-primary border border-primary-subtle", style = "font-size: 0.80rem;",
             icon("database", class = "me-1"), label),
        mod_badge
      )
    })
  })
  
  # Subpanel 2 (Distribution QC) Selectors
  # Standardized Per-Subtab Dataset Selectors
  lapply(c("sum", "scale", "box", "heat", "pca", "mds"), function(tab) {
    output[[paste0("ui_", tab, "_dataset_selector")]] <- renderUI({
      dsets <- available_bubble_datasets()
      if (is.null(dsets) || length(dsets) == 0) return(NULL)
      choices <- list()
      if (length(dsets) > 1 && tab %in% c("sum", "box", "heat", "pca", "mds")) {
        choices[[sprintf("All Datasets (Average Across %d Datasets)", length(dsets))]] <- "all_avg"
      }
      for (d in dsets) choices[[d]] <- d
      current_sel <- input[[paste0("sel_", tab, "_dataset")]] %||% choices[[1]]
      if (!isTRUE(current_sel %in% unlist(choices))) current_sel <- choices[[1]]
      selectInput(paste0("sel_", tab, "_dataset"), "Select Dataset:", choices = choices, selected = current_sel)
    })
  })

  output$ui_box_dataset_selector_cat <- renderUI({
    dsets <- available_bubble_datasets()
    if (is.null(dsets) || length(dsets) == 0) return(NULL)
    choices <- list()
    if (length(dsets) > 1) {
      choices[[sprintf("All Datasets (Average Across %d Datasets)", length(dsets))]] <- "all_avg"
    }
    for (d in dsets) choices[[d]] <- d
    current_sel <- input$sel_box_dataset %||% choices[[1]]
    if (!isTRUE(current_sel %in% unlist(choices))) current_sel <- choices[[1]]
    selectInput("sel_box_dataset", "Select Dataset:", choices = choices, selected = current_sel)
  })

  output$ui_dist_dataset_selector <- renderUI({
    dsets <- available_bubble_datasets()
    if (is.null(dsets) || length(dsets) == 0) return(NULL)
    current_sel <- input$sel_dist_dataset %||% (if (!identical(input$sel_viz_dataset_view, "all_avg") && !is.null(input$sel_viz_dataset_view)) input$sel_viz_dataset_view else dsets[1])
    if (length(current_sel) == 0 || is.null(current_sel) || !isTRUE(current_sel %in% dsets)) {
      current_sel <- dsets[1]
    }
    selectInput("sel_dist_dataset", "Dataset Reference:", choices = dsets, selected = current_sel)
  })
  output$ui_dist_modality_selector <- renderUI({
    req(rv$benchmark_df)
    has_multi <- "Modality" %in% colnames(rv$benchmark_df) &&
      (any(grepl("RNA", rv$benchmark_df$Modality, ignore.case = TRUE)) && any(grepl("ATAC", rv$benchmark_df$Modality, ignore.case = TRUE)))
    if (!has_multi) return(NULL)
    selectInput(
      "sel_dist_modality", "Assay Modality:",
      choices = c("scRNA-seq (Transcriptome)" = "RNA", "scATAC-seq (Chromatin Accessibility)" = "ATAC"),
      selected = input$sel_dist_modality %||% "RNA"
    )
  })
  output$ui_dist_sim_selector <- renderUI({
    cur_dset <- input$sel_dist_dataset %||% (input$sel_viz_dataset_view %||% "Dataset 1")
    sim_choices <- if (!is.null(rv$benchmark_df) && "Dataset" %in% colnames(rv$benchmark_df)) {
      sub_bm <- rv$benchmark_df[rv$benchmark_df$Dataset == cur_dset, ]
      if (nrow(sub_bm) > 0) unique(sub_bm$Method) else rv$methods
    } else {
      rv$methods
    }
    if (is.null(sim_choices) || length(sim_choices) == 0) sim_choices <- paste0("Simulator ", 1:6)
    current_sel <- input$sel_dist_sim %||% sim_choices[1]
    if (!current_sel %in% sim_choices) current_sel <- sim_choices[1]
    selectInput("sel_dist_sim", "Simulator:", choices = sim_choices, selected = current_sel)
  })
  
  # Standardized Per-Subtab Modality Selectors (Tabs 4, 5, 6, 7)
  lapply(c("box", "heat", "pca", "mds"), function(tab) {
    output[[paste0("ui_", tab, "_modality_selector")]] <- renderUI({
      req(rv$benchmark_df)
      choices <- c("All Modalities Combined" = "all", "scRNA-seq (Transcriptome)" = "RNA", "scATAC-seq (Chromatin Accessibility)" = "ATAC")
      current_sel <- input[[paste0("sel_", tab, "_modality")]] %||% "all"
      selectInput(paste0("sel_", tab, "_modality"), "Assay Modality:", choices = choices, selected = current_sel)
    })
  })
  output$ui_box_modality_selector_cat <- renderUI({
    req(rv$benchmark_df)
    choices <- c("All Modalities Combined" = "all", "scRNA-seq (Transcriptome)" = "RNA", "scATAC-seq (Chromatin Accessibility)" = "ATAC")
    current_sel <- input$sel_box_modality %||% "all"
    selectInput("sel_box_modality", "Assay Modality:", choices = choices, selected = current_sel)
  })
  output$ui_scale_modality_selector <- renderUI({
    req(rv$benchmark_df)
    choices <- c(
      "Combined / Total Scalability" = "total",
      "scRNA-seq Simulation Layer"  = "rna",
      "scATAC-seq Simulation Layer" = "atac"
    )
    current_sel <- input$sel_scale_modality %||% "total"
    selectInput("sel_scale_modality", "Scalability Layer:", choices = choices, selected = current_sel)
  })

  # Subpanel 8 (Cell Embeddings) Dataset & Modality Selectors
  output$ui_emb_dataset_selector <- renderUI({
    dsets <- available_bubble_datasets()
    if (is.null(dsets) || length(dsets) == 0) return(NULL)
    current_sel <- input$sel_emb_dataset %||% (if (!identical(input$sel_viz_dataset_view, "all_avg") && !is.null(input$sel_viz_dataset_view)) input$sel_viz_dataset_view else dsets[1])
    if (length(current_sel) == 0 || is.null(current_sel) || !isTRUE(current_sel %in% dsets)) {
      current_sel <- dsets[1]
    }
    selectInput("sel_emb_dataset", "Dataset to Project:", choices = dsets, selected = current_sel)
  })
  output$ui_emb_modality_selector <- renderUI({
    choices <- c("scRNA-seq (Transcriptome)" = "RNA", "scATAC-seq (Chromatin Accessibility)" = "ATAC")
    current_sel <- input$sel_emb_modality %||% "RNA"
    selectInput("sel_emb_modality", "Assay Modality:", choices = choices, selected = current_sel)
  })
  
  # ============================================================================
  # Standalone R Script Generators for Diagnostic Visualizations (Analytical Workflow)
  # ============================================================================
  generate_viz_plot_script <- function(plot_type) {
    req(rv$benchmark_df)
    df_raw <- rv$benchmark_df
    date_str <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    d_view <- input$sel_viz_dataset_view %||% (input$sel_bubble_dataset_view %||% "all_avg")
    mod_view <- input$sel_viz_modality %||% (input$sel_bubble_modality %||% "all")
    
    script_title <- switch(plot_type,
      summary = "Evaluation Summary (Radar & Horizontal Category Bars)",
      scale   = "Computational Scalability Benchmark (Runtime vs Memory)",
      box     = if (identical(input$opt_box_view_mode, "individual")) "Individual Metric Performance Barplot" else "Category Metric Distribution Boxplots",
      heat    = "Method-by-Metric Comprehensive Heatmap",
      pca     = "Principal Component Analysis (PCA) Ordination",
      mds     = "Multidimensional Scaling (MDS) Metric Space Projection",
      toupper(plot_type)
    )
    
    out <- c(
      "#!/usr/bin/env Rscript",
      "# ==============================================================================",
      sprintf("# scSimEval: %s Analysis Workflow", script_title),
      sprintf("# Generated on: %s", date_str),
      sprintf("# Active Dataset Scope: %s | Modality Scope: %s", d_view, mod_view),
      "# ==============================================================================",
      "# Standalone analytical workflow script implementing scSimEval diagnostic visualization.",
      "# Non-coding users can copy this script to HPC, cloud computing, or computing clusters",
      "# and run it directly by modifying the FILE PATH CONFIGURATION at the end of the script.",
      "# ==============================================================================",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW FUNCTION DEFINITION",
      "# ------------------------------------------------------------------------------",
      "run_visualization_workflow <- function(",
      "  benchmark_file = BENCHMARK_RESULTS_FILE,",
      "  output_dir     = OUTPUT_DIR,",
      sprintf("  dataset_scope  = '%s',", d_view),
      sprintf("  modality_scope = '%s'", mod_view),
      ") {",
      "  suppressPackageStartupMessages({",
      "    library(scSimEval)",
      "    library(ggplot2)",
      "  })",
      "",
      "  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)",
      "  cat('[scSimEval Workflow] Step 1: Environment initialized successfully.\n')",
      "",
      "  # Step 2: Ingest Benchmark Evaluation Data",
      "  if (file.exists(benchmark_file)) {",
      "    cat(sprintf('[scSimEval Workflow] Step 2: Ingesting benchmark evaluation data from: %s\n', benchmark_file))",
      "    if (endsWith(tolower(benchmark_file), '.rds')) {",
      "      obj <- readRDS(benchmark_file)",
      "      benchmark_data <- if (is.data.frame(obj)) obj else (obj$benchmark_summary_table %||% obj$benchmark_df)",
      "    } else {",
      "      benchmark_data <- utils::read.csv(benchmark_file, stringsAsFactors = FALSE, check.names = FALSE)",
      "    }",
      "  } else {",
      "    cat(sprintf('[scSimEval Workflow] Step 2: Specified file \"%s\" not found on disk.\n', benchmark_file))",
      "    cat('  Please verify or update the path in FILE PATH CONFIGURATION at the end of this script.\n')",
      "    cat('  Loading package demo benchmark database for execution verification...\n')",
      "    demo_pth <- system.file('shiny/scSimEvalApp/data/demo_benchmark_data.rds', package = 'scSimEval')",
      "    if (file.exists(demo_pth)) {",
      "      demo_obj <- readRDS(demo_pth)",
      "      benchmark_data <- demo_obj$benchmark_summary_table %||% demo_obj$benchmark_df",
      "    } else {",
      "      stop(sprintf('File \"%s\" not found. Please edit FILE PATH CONFIGURATION at the end of the script.', benchmark_file))",
      "    }",
      "  }",
      "",
      "  df <- benchmark_data",
      "  if (!'Dataset' %in% colnames(df)) df$Dataset <- 'Dataset 1'",
      "  if (exists('standardize_benchmark_categories', where = asNamespace('scSimEval'), mode = 'function')) {",
      "    df <- scSimEval::standardize_benchmark_categories(df)",
      "  }",
      "",
      "  # Step 3: Modality & Multi-Dataset Harmonization",
      "  if ('Modality' %in% colnames(df) && !identical(modality_scope, 'all')) {",
      "    cat(sprintf('[scSimEval Workflow] Filtering for modality: %s\n', modality_scope))",
      "    if (identical(modality_scope, 'RNA')) df <- df[df$Modality %in% c('RNA', 'Computation'), , drop = FALSE]",
      "    else if (identical(modality_scope, 'ATAC')) df <- df[df$Modality %in% c('ATAC', 'Computation'), , drop = FALSE]",
      "    else if (identical(modality_scope, 'Joint')) df <- df[df$Modality %in% c('Joint (Cross-Modal)', 'Computation'), , drop = FALSE]",
      "  }",
      "",
      "  all_dsets <- unique(as.character(df$Dataset))",
      "  all_dsets <- all_dsets[!is.na(all_dsets) & all_dsets != '']",
      "  if (identical(dataset_scope, 'all_avg') && length(all_dsets) > 1) {",
      "    score_col <- if ('Score' %in% colnames(df)) 'Score' else 'Value'",
      "    cat(sprintf('[scSimEval Workflow] Computing cross-dataset average across %d datasets...\n', length(all_dsets)))",
      "    agg_score <- stats::aggregate(df[[score_col]], by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)",
      "    colnames(agg_score)[4] <- score_col",
      "    agg_score$Dataset <- sprintf('All Datasets (Average of %d Datasets)', length(all_dsets))",
      "    if ('Score_Raw' %in% colnames(df)) agg_score$Score_Raw <- stats::aggregate(df$Score_Raw, by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)[[4]]",
      "    if ('Score_Norm' %in% colnames(df)) agg_score$Score_Norm <- stats::aggregate(df$Score_Norm, by = list(Method = df$Method, Category = df$Category, Metric = df$Metric), FUN = mean, na.rm = TRUE)[[4]]",
      "    df <- agg_score",
      "  } else if (!identical(dataset_scope, 'all_avg') && dataset_scope %in% all_dsets) {",
      "    df <- df[df$Dataset == dataset_scope, , drop = FALSE]",
      "  }",
      ""
    )
    
    out_body <- switch(plot_type,
      "summary" = c(
        "  # Step 4: Leaderboard Ranking & Evaluation Summary Figure",
        "  cat_summary <- stats::aggregate(df$Score, by = list(Method = df$Method, Category = df$Category), FUN = mean, na.rm = TRUE)",
        "  colnames(cat_summary)[3] <- 'Mean_Category_Score'",
        "  overall_ranks <- stats::aggregate(cat_summary$Mean_Category_Score, by = list(Method = cat_summary$Method), FUN = mean, na.rm = TRUE)",
        "  colnames(overall_ranks) <- c('Method', 'Overall_Composite_Score')",
        "  overall_ranks <- overall_ranks[order(-overall_ranks$Overall_Composite_Score), ]",
        "  cat('[scSimEval Workflow] Overall Simulator Ranking Summary:\n')",
        "  print(overall_ranks)",
        "",
        "  cat('[scSimEval Workflow] Step 5: Generating horizontal summary visualization...\n')",
        sprintf("  p <- plot_evaluation_summary(data = df, show_labels = %s, normalize_scores = %s, base_size = 14)", input$chk_sum_labels %||% TRUE, input$chk_sum_norm %||% TRUE),
        "",
        "  pdf_file  <- file.path(output_dir, 'scSimEval_evaluation_summary.pdf')",
        "  jpeg_file <- file.path(output_dir, 'scSimEval_evaluation_summary.jpeg')",
        "  ggplot2::ggsave(pdf_file, plot = p, width = 13, height = 7.5, device = 'pdf')",
        "  ggplot2::ggsave(jpeg_file, plot = p, width = 13, height = 7.5, dpi = 600)",
        "  cat(sprintf('[scSimEval Workflow] Figures successfully saved:\n  - Vector PDF: %s\n  - Publication JPEG: %s\n', pdf_file, jpeg_file))"
      ),
      "scale" = c(
        "  # Step 4: Computational Efficiency & Scalability Dashboard",
        "  cat('[scSimEval Workflow] Generating scalability benchmark dashboard...\n')",
        sprintf("  p <- plot_scalability_benchmark(benchmark_data = df, type = '%s', base_size = 13)", input$sel_scale_type %||% "composite"),
        "",
        "  pdf_file  <- file.path(output_dir, 'scSimEval_scalability_benchmark.pdf')",
        "  jpeg_file <- file.path(output_dir, 'scSimEval_scalability_benchmark.jpeg')",
        "  ggplot2::ggsave(pdf_file, plot = p, width = 13, height = 8, device = 'pdf')",
        "  ggplot2::ggsave(jpeg_file, plot = p, width = 13, height = 8, dpi = 600)",
        "  cat(sprintf('[scSimEval Workflow] Figures successfully saved:\n  - Vector PDF: %s\n  - Publication JPEG: %s\n', pdf_file, jpeg_file))"
      ),
      "box" = if (identical(input$opt_box_view_mode, "individual")) {
        c(
          sprintf("  metric_name <- '%s'", input$sel_box_metric_single %||% "Mean Expression"),
          sprintf("  score_scale <- '%s'", input$sel_box_indiv_score %||% "normalized"),
          "  cat(sprintf('[scSimEval Workflow] Step 4: Generating ranked barplot for: %s (%s)...\n', metric_name, score_scale))",
          "  p <- plot_individual_metric_bar(benchmark_data = df, metric = metric_name, score_type = score_scale, base_size = 12)",
          "",
          "  pdf_file  <- file.path(output_dir, 'scSimEval_metric_bar.pdf')",
          "  jpeg_file <- file.path(output_dir, 'scSimEval_metric_bar.jpeg')",
          "  ggplot2::ggsave(pdf_file, plot = p, width = 9.5, height = 6.5, device = 'pdf')",
          "  ggplot2::ggsave(jpeg_file, plot = p, width = 9.5, height = 6.5, dpi = 600)",
          "  cat(sprintf('[scSimEval Workflow] Figures successfully saved:\n  - Vector PDF: %s\n  - Publication JPEG: %s\n', pdf_file, jpeg_file))"
        )
      } else {
        c(
          sprintf("  cat_subset  <- %s", if (identical(input$sel_box_cat_group, "all")) "NULL" else paste0("'", input$sel_box_cat_group, "'")),
          sprintf("  score_scale <- '%s'", input$sel_box_score_type %||% "normalized"),
          "  cat('[scSimEval Workflow] Step 4: Generating category distribution boxplots with jittered points...\n')",
          "  p <- plot_metric_boxplots(benchmark_data = df, categories = cat_subset, score_type = score_scale, facet_by = 'category', base_size = 11)",
          "",
          "  pdf_file  <- file.path(output_dir, 'scSimEval_metric_boxplots.pdf')",
          "  jpeg_file <- file.path(output_dir, 'scSimEval_metric_boxplots.jpeg')",
          "  ggplot2::ggsave(pdf_file, plot = p, width = 12, height = 8.5, device = 'pdf')",
          "  ggplot2::ggsave(jpeg_file, plot = p, width = 12, height = 8.5, dpi = 600)",
          "  cat(sprintf('[scSimEval Workflow] Figures successfully saved:\n  - Vector PDF: %s\n  - Publication JPEG: %s\n', pdf_file, jpeg_file))"
        )
      },
      "heat" = c(
        sprintf("  cat_filter  <- %s", if (identical(input$sel_heat_cat, "all")) "NULL" else paste0("'", input$sel_heat_cat, "'")),
        sprintf("  scale_type  <- '%s'", if (!is.null(input$sel_heat_score_type) && identical(input$sel_heat_score_type, "raw")) "raw" else "relative"),
        "  cat('[scSimEval Workflow] Step 4: Generating publication metric heatmap...\n')",
        "  p <- plot_metric_heatmap(benchmark_data = df, category = cat_filter, scale_fill = scale_type, cluster_rows = FALSE, cluster_cols = FALSE, base_size = 11)",
        "",
        sprintf("  w_in <- %.2f", (input$sld_heat_width %||% 880) / 90),
        sprintf("  h_in <- %.2f", (input$sld_heat_height %||% 1100) / 90),
        "  pdf_file  <- file.path(output_dir, 'scSimEval_metric_heatmap.pdf')",
        "  jpeg_file <- file.path(output_dir, 'scSimEval_metric_heatmap.jpeg')",
        "  ggplot2::ggsave(pdf_file, plot = p, width = w_in, height = h_in, device = 'pdf')",
        "  ggplot2::ggsave(jpeg_file, plot = p, width = w_in, height = h_in, dpi = 600)",
        "  cat(sprintf('[scSimEval Workflow] Figures successfully saved:\n  - Vector PDF: %s\n  - Publication JPEG: %s\n', pdf_file, jpeg_file))"
      ),
      "pca" = c(
        sprintf("  cat_filter   <- %s", if (identical(input$sel_pca_cat, "all")) "NULL" else paste0("'", input$sel_pca_cat, "'")),
        sprintf("  panel_type   <- '%s'", input$sel_pca_panel %||% "both"),
        sprintf("  top_loadings <- %d", input$num_pca_top_metrics %||% 8),
        "  cat('[scSimEval Workflow] Step 4: Generating PCA ordination and loadings plot...\n')",
        "  p <- plot_metric_pca(benchmark_data = df, category = cat_filter, panel = panel_type, top_n_loadings = top_loadings, base_size = 12)",
        "",
        sprintf("  w_in <- %.2f", (input$sld_pca_width %||% 950) / 90),
        sprintf("  h_in <- %.2f", (input$sld_pca_height %||% 1100) / 90),
        "  pdf_file  <- file.path(output_dir, 'scSimEval_metric_pca.pdf')",
        "  jpeg_file <- file.path(output_dir, 'scSimEval_metric_pca.jpeg')",
        "  ggplot2::ggsave(pdf_file, plot = p, width = w_in, height = h_in, device = 'pdf')",
        "  ggplot2::ggsave(jpeg_file, plot = p, width = w_in, height = h_in, dpi = 600)",
        "  cat(sprintf('[scSimEval Workflow] Figures successfully saved:\n  - Vector PDF: %s\n  - Publication JPEG: %s\n', pdf_file, jpeg_file))"
      ),
      "mds" = c(
        sprintf("  cat_filter <- %s", if (identical(input$sel_mds_cat, "all")) "NULL" else paste0("'", input$sel_mds_cat, "'")),
        sprintf("  ord_by     <- '%s'", input$sel_mds_by %||% "simulators"),
        "  cat('[scSimEval Workflow] Step 4: Generating MDS metric space projection...\n')",
        "  p <- plot_metric_mds(benchmark_data = df, category = cat_filter, ordination_by = ord_by, base_size = 13)",
        "",
        sprintf("  w_in <- %.2f", (input$sld_mds_width %||% 900) / 80),
        sprintf("  h_in <- %.2f", (input$sld_mds_height %||% 620) / 80),
        "  pdf_file  <- file.path(output_dir, 'scSimEval_metric_mds.pdf')",
        "  jpeg_file <- file.path(output_dir, 'scSimEval_metric_mds.jpeg')",
        "  ggplot2::ggsave(pdf_file, plot = p, width = w_in, height = h_in, device = 'pdf')",
        "  ggplot2::ggsave(jpeg_file, plot = p, width = w_in, height = h_in, dpi = 600)",
        "  cat(sprintf('[scSimEval Workflow] Figures successfully saved:\n  - Vector PDF: %s\n  - Publication JPEG: %s\n', pdf_file, jpeg_file))"
      )
    )
    
    out <- c(out,
      out_body,
      "}",
      "",
      get_workflow_file_paths_block(plot_type),
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW INVOCATION",
      "# ------------------------------------------------------------------------------",
      sprintf("run_visualization_workflow(benchmark_file = BENCHMARK_RESULTS_FILE, output_dir = OUTPUT_DIR, dataset_scope = '%s', modality_scope = '%s')", d_view, mod_view),
      ""
    )
    
    return(paste(out, collapse = "\n"))
  }
  
  # Standalone Distribution QC Analysis Workflow Script Generator
  generate_dist_qc_script <- function() {
    cur_dset <- input$sel_dist_dataset %||% (input$sel_viz_dataset_view %||% "Dataset 1")
    if (identical(cur_dset, "all_avg")) {
      dsets <- available_bubble_datasets()
      cur_dset <- if (length(dsets) > 0) dsets[1] else "Dataset 1"
    }
    sim_choice <- input$sel_dist_sim %||% "Simulator 1"
    mod_choice <- input$sel_dist_modality %||% "RNA"
    date_str <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    lay <- input$sel_dist_layout %||% "full"
    
    out <- c(
      "#!/usr/bin/env Rscript",
      "# ==============================================================================",
      "# scSimEval: Comparative Distribution QC Analysis Workflow",
      sprintf("# Generated on: %s", date_str),
      sprintf("# Reference Dataset: %s | Evaluated Simulator: %s | Modality: %s", cur_dset, sim_choice, mod_choice),
      "# ==============================================================================",
      "# Standalone analytical workflow for evaluating 14 statistical distribution properties",
      "# comparing empirical reference counts against simulated single-cell profiles.",
      "# Non-coding users can copy this script to HPC, cloud computing, or computing clusters",
      "# and run it directly by modifying the FILE PATH CONFIGURATION at the end of the script.",
      "# ==============================================================================",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW STEP 1: Environment Setup & Matrix Loader Helper",
      "# ------------------------------------------------------------------------------",
      "suppressPackageStartupMessages({",
      "  library(scSimEval)",
      "  library(Matrix)",
      "  library(ggplot2)",
      "})",
      "",
      "load_matrix_safe <- function(path, placeholder = 'matrix') {",
      "  if (is.null(path) || is.na(path) || path == '') stop(sprintf('Path for %s is empty.', placeholder))",
      "  if (!file.exists(path)) {",
      "    warning(sprintf('File \"%s\" not found. Using small placeholder matrix for pipeline verification.', path))",
      "    set.seed(42)",
      "    m <- Matrix::rsparsematrix(nrow = 100, ncol = 50, density = 0.2)",
      "    rownames(m) <- paste0('Feature_', seq_len(100))",
      "    colnames(m) <- paste0('Cell_', seq_len(50))",
      "    return(abs(m))",
      "  }",
      "  if (endsWith(tolower(path), '.rds')) {",
      "    m <- readRDS(path)",
      "    if (inherits(m, 'SingleCellExperiment')) m <- SummarizedExperiment::assay(m, 'counts')",
      "    else if (inherits(m, 'Seurat')) m <- tryCatch(m[['RNA']]$counts, error = function(e) tryCatch(Seurat::GetAssayData(m, slot = 'counts'), error = function(e2) m))",
      "    else if (is.list(m) && !is.matrix(m) && !inherits(m, 'Matrix')) {",
      "      for (item in m) if (is.matrix(item) || inherits(item, 'Matrix')) { m <- item; break }",
      "    }",
      "  } else {",
      "    sep_char <- if (endsWith(tolower(path), '.csv')) ',' else '\t'",
      "    df <- utils::read.table(path, header = TRUE, row.names = 1, sep = sep_char, check.names = FALSE)",
      "    m <- as.matrix(df)",
      "  }",
      "  if (!inherits(m, 'dgCMatrix')) m <- Matrix::Matrix(m, sparse = TRUE)",
      "  return(m)",
      "}",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW FUNCTION DEFINITION",
      "# ------------------------------------------------------------------------------",
      "run_dist_qc_workflow <- function(",
      "  ref_file   = REF_DATA_FILE,",
      "  sim_file   = SIM_DATA_FILE,",
      sprintf("  layout     = '%s',", lay),
      "  output_dir = OUTPUT_DIR",
      ") {",
      "  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)",
      "  cat('[scSimEval Workflow] Step 2: Loading reference and simulated count matrices...\n')",
      "  ref_data <- load_matrix_safe(ref_file, placeholder = 'Reference Matrix')",
      "  sim_data <- load_matrix_safe(sim_file, placeholder = 'Simulated Matrix')",
      "",
      "  cat(sprintf('[scSimEval Workflow] Step 3: Reference: %d features x %d cells | Simulated: %d features x %d cells\n', nrow(ref_data), ncol(ref_data), nrow(sim_data), ncol(sim_data)))",
      "  common_features <- intersect(rownames(ref_data), rownames(sim_data))",
      "  if (length(common_features) > 0 && length(common_features) < nrow(ref_data)) {",
      "    ref_data <- ref_data[common_features, , drop = FALSE]",
      "    sim_data <- sim_data[common_features, , drop = FALSE]",
      "  }",
      "",
      "  cat('[scSimEval Workflow] Step 4: Generating 14-panel Distribution QC diagnostic figure...\n')",
      "  p <- plot_distribution_qc(ref_data = ref_data, sim_data = sim_data, layout = layout, base_size = 13)",
      "",
      "  pdf_file  <- file.path(output_dir, 'scSimEval_distribution_qc.pdf')",
      "  jpeg_file <- file.path(output_dir, 'scSimEval_distribution_qc.jpeg')",
      "  ggplot2::ggsave(pdf_file, plot = p, width = 14, height = 9, device = 'pdf')",
      "  ggplot2::ggsave(jpeg_file, plot = p, width = 14, height = 9, dpi = 600)",
      "  cat(sprintf('[scSimEval Workflow] Figures successfully saved:\n  - Vector PDF: %s\n  - Publication JPEG: %s\n', pdf_file, jpeg_file))",
      "}",
      "",
      get_workflow_file_paths_block("dist_qc"),
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW INVOCATION",
      "# ------------------------------------------------------------------------------",
      "run_dist_qc_workflow(ref_file = REF_DATA_FILE, sim_file = SIM_DATA_FILE, output_dir = OUTPUT_DIR)",
      ""
    )
    return(paste(out, collapse = "\n"))
  }
  
  # Standalone Cell Embeddings Analysis Workflow Script Generator
  generate_cell_emb_script <- function() {
    cur_dset <- input$sel_emb_dataset %||% (input$sel_viz_dataset_view %||% "Dataset 1")
    if (identical(cur_dset, "all_avg")) {
      dsets <- available_bubble_datasets()
      cur_dset <- if (length(dsets) > 0) dsets[1] else "Dataset 1"
    }
    date_str <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    red <- input$sel_emb_reduction %||% "umap"
    lay <- input$sel_emb_layout %||% "facet"
    col_by <- input$sel_emb_color %||% "cell_type"
    lay_arg <- if (identical(lay, "side_by_side") || identical(lay, "compare")) "side_by_side" else "facet"
    
    out <- c(
      "#!/usr/bin/env Rscript",
      "# ==============================================================================",
      "# scSimEval: Cell Embeddings & Latent Space Concordance Workflow",
      sprintf("# Generated on: %s", date_str),
      sprintf("# Dataset: %s | Reduction: %s | Layout: %s", cur_dset, red, lay_arg),
      "# ==============================================================================",
      "# Standalone analytical workflow for projecting reference and simulated single cells into a",
      "# unified latent coordinate space (UMAP/t-SNE/PCA) and quantifying cluster fidelity.",
      "# Non-coding users can copy this script to HPC, cloud computing, or computing clusters",
      "# and run it directly by modifying the FILE PATH CONFIGURATION at the end of the script.",
      "# ==============================================================================",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW STEP 1: Environment Setup & Helpers",
      "# ------------------------------------------------------------------------------",
      "suppressPackageStartupMessages({",
      "  library(scSimEval)",
      "  library(Matrix)",
      "  library(ggplot2)",
      "})",
      "",
      "load_matrix_safe <- function(path, placeholder = 'matrix') {",
      "  if (is.null(path) || is.na(path) || path == '') stop(sprintf('Path for %s is empty.', placeholder))",
      "  if (!file.exists(path)) {",
      "    warning(sprintf('File \"%s\" not found. Using small placeholder matrix for pipeline verification.', path))",
      "    set.seed(42)",
      "    m <- Matrix::rsparsematrix(nrow = 100, ncol = 50, density = 0.2)",
      "    rownames(m) <- paste0('Feature_', seq_len(100))",
      "    colnames(m) <- paste0('Cell_', seq_len(50))",
      "    return(abs(m))",
      "  }",
      "  if (endsWith(tolower(path), '.rds')) {",
      "    m <- readRDS(path)",
      "    if (inherits(m, 'SingleCellExperiment')) m <- SummarizedExperiment::assay(m, 'counts')",
      "    else if (inherits(m, 'Seurat')) m <- tryCatch(m[['RNA']]$counts, error = function(e) tryCatch(Seurat::GetAssayData(m, slot = 'counts'), error = function(e2) m))",
      "    else if (is.list(m) && !is.matrix(m) && !inherits(m, 'Matrix')) {",
      "      for (item in m) if (is.matrix(item) || inherits(item, 'Matrix')) { m <- item; break }",
      "    }",
      "  } else {",
      "    sep_char <- if (endsWith(tolower(path), '.csv')) ',' else '\t'",
      "    df <- utils::read.table(path, header = TRUE, row.names = 1, sep = sep_char, check.names = FALSE)",
      "    m <- as.matrix(df)",
      "  }",
      "  if (!inherits(m, 'dgCMatrix')) m <- Matrix::Matrix(m, sparse = TRUE)",
      "  return(m)",
      "}",
      "",
      "load_vector_safe <- function(path) {",
      "  if (is.null(path) || is.na(path) || path == '' || !file.exists(path)) return(NULL)",
      "  lower_p <- tolower(path)",
      "  if (endsWith(lower_p, '.rds')) {",
      "    v <- readRDS(path)",
      "    if (is.data.frame(v) || is.matrix(v)) v <- v[, 1]",
      "    return(as.character(v))",
      "  }",
      "  if (endsWith(lower_p, '.csv') || endsWith(lower_p, '.tsv')) {",
      "    sep_char <- if (endsWith(lower_p, '.csv')) ',' else '\t'",
      "    df <- tryCatch(utils::read.table(path, header = TRUE, sep = sep_char, stringsAsFactors = FALSE), error = function(e) NULL)",
      "    if (!is.null(df) && is.data.frame(df) && ncol(df) >= 1) return(as.character(df[[ncol(df)]]))",
      "  }",
      "  return(as.character(readLines(path)))",
      "}",
      "",
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW FUNCTION DEFINITION",
      "# ------------------------------------------------------------------------------",
      "run_cell_embeddings_workflow <- function(",
      "  ref_file        = REF_MATRIX_FILE,",
      "  sim_files       = SIMULATOR_FILES,",
      "  cell_types_file = CELL_TYPES_FILE,",
      "  batch_file      = BATCH_FILE,",
      sprintf("  reduction       = '%s',", red),
      sprintf("  layout          = '%s',", lay_arg),
      sprintf("  color_by        = '%s',", col_by),
      "  output_dir      = OUTPUT_DIR",
      ") {",
      "  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)",
      "  cat('[scSimEval Workflow] Step 2: Loading reference and simulated count matrices...\n')",
      "  ref_mat <- load_matrix_safe(ref_file, placeholder = 'Reference Matrix')",
      "  sim_mats <- list()",
      "  if (is.list(sim_files)) {",
      "    for (nm in names(sim_files)) {",
      "      f <- sim_files[[nm]]",
      "      if (is.list(f)) f <- f$rna %||% f[[1]]",
      "      sim_mats[[nm]] <- load_matrix_safe(f, placeholder = paste0('Simulated Matrix: ', nm))",
      "    }",
      "  }",
      "  cell_types <- load_vector_safe(cell_types_file)",
      "  batch      <- load_vector_safe(batch_file)",
      "",
      "  cat(sprintf('[scSimEval Workflow] Step 3: Computing unified %s cell embeddings...\n', toupper(reduction)))",
      "  emb_obj <- compute_dataset_embeddings(",
      "    reference  = ref_mat,",
      "    simulated  = sim_mats,",
      "    reduction  = reduction,",
      "    cell_types = cell_types,",
      "    batch      = batch",
      "  )",
      "",
      "  cat('[scSimEval Workflow] Step 4: Quantitative Embedding Quality Metrics:\n')",
      "  quality_metrics <- compute_embedding_quality_metrics(emb_obj)",
      "  print(quality_metrics)",
      "",
      "  cat('[scSimEval Workflow] Step 5: Rendering cell embeddings visualization...\n')",
      "  p <- plot_dataset_embeddings(embedding_data = emb_obj, reduction = reduction, layout = layout, color_by = color_by)",
      "",
      "  pdf_file  <- file.path(output_dir, 'scSimEval_cell_embeddings.pdf')",
      "  jpeg_file <- file.path(output_dir, 'scSimEval_cell_embeddings.jpeg')",
      "  ggplot2::ggsave(pdf_file, plot = p, width = 13, height = 8.5, device = 'pdf')",
      "  ggplot2::ggsave(jpeg_file, plot = p, width = 13, height = 8.5, dpi = 600)",
      "  cat(sprintf('[scSimEval Workflow] Figures successfully saved:\n  - Vector PDF: %s\n  - Publication JPEG: %s\n', pdf_file, jpeg_file))",
      "}",
      "",
      get_workflow_file_paths_block("cell_emb"),
      "# ------------------------------------------------------------------------------",
      "# WORKFLOW INVOCATION",
      "# ------------------------------------------------------------------------------",
      "run_cell_embeddings_workflow(ref_file = REF_MATRIX_FILE, sim_files = SIMULATOR_FILES, cell_types_file = CELL_TYPES_FILE, batch_file = BATCH_FILE, output_dir = OUTPUT_DIR)",
      ""
    )
    return(paste(out, collapse = "\n"))
  }

  # Download Handlers for R Scripts across all 8 sub-panels
  output$download_sum_r_script <- downloadHandler(
    filename = function() { paste0("scSimEval_evaluation_summary_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("summary"), file) }
  )
  output$download_sum_r_script_header <- downloadHandler(
    filename = function() { paste0("scSimEval_evaluation_summary_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("summary"), file) }
  )
  output$download_dist_r_script <- downloadHandler(
    filename = function() { paste0("scSimEval_distribution_qc_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_dist_qc_script(), file) }
  )
  output$download_dist_r_script_header <- downloadHandler(
    filename = function() { paste0("scSimEval_distribution_qc_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_dist_qc_script(), file) }
  )
  output$download_scale_r_script <- downloadHandler(
    filename = function() { paste0("scSimEval_scalability_benchmark_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("scale"), file) }
  )
  output$download_scale_r_script_header <- downloadHandler(
    filename = function() { paste0("scSimEval_scalability_benchmark_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("scale"), file) }
  )
  output$download_box_r_script <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_plot_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("box"), file) }
  )
  output$download_box_r_script_header <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_plot_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("box"), file) }
  )
  output$download_box_cat_r_script <- downloadHandler(
    filename = function() { paste0("scSimEval_category_metric_plot_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("box"), file) }
  )
  output$download_heat_r_script <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_heatmap_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("heat"), file) }
  )
  output$download_heat_r_script_header <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_heatmap_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("heat"), file) }
  )
  output$download_pca_r_script <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_pca_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("pca"), file) }
  )
  output$download_pca_r_script_header <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_pca_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("pca"), file) }
  )
  output$download_mds_r_script <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_mds_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("mds"), file) }
  )
  output$download_mds_r_script_header <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_mds_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_viz_plot_script("mds"), file) }
  )
  output$download_emb_r_script <- downloadHandler(
    filename = function() { paste0("scSimEval_cell_embeddings_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_cell_emb_script(), file) }
  )
  output$download_emb_r_script_header <- downloadHandler(
    filename = function() { paste0("scSimEval_cell_embeddings_", Sys.Date(), ".R") },
    content = function(file) { writeLines(generate_cell_emb_script(), file) }
  )
  
  # 1. Evaluation Summary
  output$ui_plot_eval_summary <- renderUI({
    if (rv_plot_triggers$sum == 0) {
      return(render_plot_placeholder(
        title = "Evaluation Summary Ready to Plot",
        desc = "Click below to confirm and render the radar and category-averaged performance summary bars.",
        button_id = "btn_plot_sum_inline",
        btn_label = "Confirm & Plot Evaluation Summary",
        icon_name = "chart-pie"
      ))
    }
    plotOutput("plot_eval_summary", height = "540px")
  })
  
  eval_summary_reactive <- reactive({
    req(rv$benchmark_df)
    req(rv_plot_triggers$sum > 0)
    df <- get_viz_benchmark_data(dataset_choice = input$sel_sum_dataset)
    withProgress(message = "Rendering Evaluation Summary", detail = "Aggregating category scores...", {
      plot_evaluation_summary(
        data = df,
        show_labels = input$chk_sum_labels,
        normalize_scores = input$chk_sum_norm,
        base_size = 14
      )
    })
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
    if (rv_plot_triggers$dist == 0) {
      return(render_plot_placeholder(
        title = "Distribution QC Comparison Ready",
        desc = "Click below to confirm and compare gene expression, library size, and dropout density curves between reference and simulated datasets.",
        button_id = "btn_plot_dist_inline",
        btn_label = "Confirm & Plot Distribution QC",
        icon_name = "chart-line"
      ))
    }
    plot_h <- if (identical(input$sel_dist_layout, "density_only") || identical(input$sel_dist_layout, "density")) "550px" else "850px"
    plotOutput("plot_dist_qc", height = plot_h)
  })
  dist_qc_reactive <- reactive({
    req(rv_plot_triggers$dist > 0)
    cur_dset <- input$sel_dist_dataset %||% (input$sel_viz_dataset_view %||% "Dataset 1")
    if (identical(cur_dset, "all_avg")) cur_dset <- "Dataset 1"
    cur_mod <- input$sel_dist_modality %||% "RNA"
    cur_dset_clean <- sub(" \\(.*\\)$", "", cur_dset)
    
    # Try resolving reference matrix
    ref_mat <- NULL
    ref_try_keys <- c(
      paste0(cur_dset, "_", cur_mod),
      paste0(cur_dset_clean, "_", cur_mod),
      paste0(cur_dset, " (sc", cur_mod, "-seq)"),
      paste0(cur_dset_clean, " (sc", cur_mod, "-seq)"),
      cur_dset,
      cur_dset_clean,
      paste0(cur_dset, "_", tolower(cur_mod)),
      paste0(cur_dset_clean, "_", tolower(cur_mod))
    )
    for (rk in ref_try_keys) {
      if (!is.null(rv$ref_matrices) && rk %in% names(rv$ref_matrices)) { ref_mat <- rv$ref_matrices[[rk]]; break }
      if (!is.null(rv$all_matrices) && paste0("[Reference] ", rk) %in% names(rv$all_matrices)) { ref_mat <- rv$all_matrices[[paste0("[Reference] ", rk)]]; break }
    }
    if (is.null(ref_mat) && !is.null(rv$ref_matrices) && length(rv$ref_matrices) > 0) {
      mod_hits <- names(rv$ref_matrices)[grepl(cur_mod, names(rv$ref_matrices), ignore.case = TRUE)]
      if (length(mod_hits) > 0) {
        dset_hits <- mod_hits[grepl(cur_dset_clean, mod_hits, fixed = TRUE)]
        ref_mat <- if (length(dset_hits) > 0) rv$ref_matrices[[dset_hits[1]]] else rv$ref_matrices[[mod_hits[1]]]
      } else if (cur_dset %in% names(rv$ref_matrices)) {
        ref_mat <- rv$ref_matrices[[cur_dset]]
      } else {
        ref_mat <- rv$ref_matrices[[1]]
      }
    }
    if (is.null(ref_mat)) ref_mat <- rv$toy_ref
    req(ref_mat)
    
    sim_choice <- input$sel_dist_sim %||% "Simulator 1"
    sim_mat <- NULL
    try_keys <- c(
      paste0(sim_choice, "_", cur_mod, "_", cur_dset),
      paste0(sim_choice, "_", cur_mod, "_", cur_dset_clean),
      paste0(sim_choice, "_", tolower(cur_mod), "_", cur_dset),
      paste0(sim_choice, "_", tolower(cur_mod), "_", cur_dset_clean),
      paste0(sim_choice, "_", cur_dset),
      paste0(sim_choice, "_", cur_dset_clean),
      paste0(sim_choice, "_", cur_mod),
      paste0(sim_choice, " (", cur_dset, ")"),
      paste0(sim_choice, " (", cur_dset_clean, ")"),
      paste0("Simulated: ", sim_choice, " (", cur_dset, ")"),
      paste0("[Simulated] ", sim_choice, " ", cur_mod, " (", cur_dset_clean, ")"),
      paste0("[Simulated] ", sim_choice, "_", cur_dset),
      paste0("[", cur_dset, "] Simulated: ", sim_choice),
      paste0("[", cur_dset, "] Simulated ", cur_mod, ": ", sim_choice),
      sim_choice,
      paste0("Simulated: ", sim_choice)
    )
    for (k in try_keys) {
      if (!is.null(rv$sim_matrices) && k %in% names(rv$sim_matrices)) { sim_mat <- rv$sim_matrices[[k]]; break }
      if (!is.null(rv$all_matrices) && k %in% names(rv$all_matrices)) { sim_mat <- rv$all_matrices[[k]]; break }
    }
    if (is.null(sim_mat) && !is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) {
      matches <- names(rv$sim_matrices)[grepl(paste0("^", sim_choice, "(_| |\\()"), names(rv$sim_matrices))]
      if (length(matches) > 0) {
        mod_m <- matches[grepl(cur_mod, matches, ignore.case = TRUE)]
        sim_mat <- if (length(mod_m) > 0) rv$sim_matrices[[mod_m[1]]] else rv$sim_matrices[[matches[1]]]
      }
    }
    if (is.null(sim_mat) && !is.null(rv$all_matrices) && length(rv$all_matrices) > 0) {
      matches <- names(rv$all_matrices)[grepl(sim_choice, names(rv$all_matrices), fixed = TRUE)]
      if (length(matches) > 0) {
        mod_m <- matches[grepl(cur_mod, matches, ignore.case = TRUE)]
        sim_mat <- if (length(mod_m) > 0) rv$all_matrices[[mod_m[1]]] else rv$all_matrices[[matches[1]]]
      }
    }
    if (is.null(sim_mat) && !is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) {
      sim_mat <- rv$sim_matrices[[1]]
    }
    if (is.null(sim_mat)) sim_mat <- rv$toy_sim
    req(sim_mat)
    
    withProgress(message = "Rendering Distribution QC", detail = "Evaluating expression densities...", {
      lay <- if (!is.null(input$sel_dist_layout) && input$sel_dist_layout %in% c("full", "density_only", "grid_compact")) {
        input$sel_dist_layout
      } else if (identical(input$sel_dist_layout, "density")) {
        "density_only"
      } else {
        "full"
      }
      plot_distribution_qc(
        ref_data = ref_mat,
        sim_data = sim_mat,
        layout   = lay,
        base_size = 13
      )
    })
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
  output$ui_plot_scale_bench <- renderUI({
    if (rv_plot_triggers$scale == 0) {
      return(render_plot_placeholder(
        title = "Scalability Benchmark Ready to Plot",
        desc = "Click below to confirm and visualize runtime (seconds) and peak memory usage (MiB) across simulators.",
        button_id = "btn_plot_scale_inline",
        btn_label = "Confirm & Plot Scalability",
        icon_name = "microchip"
      ))
    }
    plotOutput("plot_scale_bench", height = "620px")
  })
  scale_bench_reactive <- reactive({
    req(rv$benchmark_df)
    req(rv_plot_triggers$scale > 0)
    scale_mod <- input$sel_scale_modality %||% "total"
    df <- get_viz_benchmark_data(
      dataset_choice = input$sel_scale_dataset,
      modality_choice = if (scale_mod == "rna") "RNA" else if (scale_mod == "atac") "ATAC" else "all"
    )
    withProgress(message = "Rendering Scalability Benchmark", detail = "Extracting computational metrics...", {
      plot_scalability_benchmark(
        benchmark_data = df,
        type = input$sel_scale_type,
        modality = scale_mod,
        base_size = 13
      )
    })
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
  
  # 4. Metric Plots
  output$ui_box_metric_picker <- renderUI({
    df <- get_viz_benchmark_data(dataset_choice = input$sel_box_dataset, modality_choice = input$sel_box_modality)
    req(input$sel_box_cat_first)
    sub_df <- df[df$Category == input$sel_box_cat_first, , drop = FALSE]
    avail_metrics <- sort(unique(sub_df$Metric))
    selectInput("sel_box_metric_single", "2. Choose Metric Name:", choices = avail_metrics, selected = avail_metrics[1])
  })
  
  metric_box_reactive <- reactive({
    req(rv$benchmark_df)
    req(rv_plot_triggers$box > 0)
    df <- get_viz_benchmark_data(
      dataset_choice = if (identical(input$opt_box_view_mode, "category")) input$sel_box_dataset_cat else input$sel_box_dataset,
      modality_choice = if (identical(input$opt_box_view_mode, "category")) input$sel_box_modality_cat else input$sel_box_modality
    )
    withProgress(message = "Rendering Metric Plots", detail = "Formatting metric scores...", {
      if (identical(input$opt_box_view_mode, "individual")) {
        req(input$sel_box_metric_single)
        score_type_val <- if (!is.null(input$sel_box_indiv_score)) input$sel_box_indiv_score else "normalized"
        plot_individual_metric_bar(
          benchmark_data = df,
          metric         = input$sel_box_metric_single,
          score_type     = score_type_val,
          base_size      = 12
        )
      } else {
        if (identical(input$sel_box_cat_group, "all")) {
          plot_metric_boxplots(
            benchmark_data = df,
            categories     = NULL,
            score_type     = input$sel_box_score_type,
            facet_by       = "category",
            base_size      = 11
          )
        } else {
          plot_category_metric_bars(
            benchmark_data = df,
            category       = input$sel_box_cat_group,
            score_type     = input$sel_box_score_type,
            base_size      = 11
          )
        }
      }
    })
  })
  
  output$ui_plot_metric_boxes <- renderUI({
    if (rv_plot_triggers$box == 0) {
      return(render_plot_placeholder(
        title = "Metric Plots Ready to Generate",
        desc = "Choose your metric or category group, then click below to confirm and render comparison plots.",
        button_id = "btn_plot_box_inline",
        btn_label = "Confirm & Plot Metric",
        icon_name = "chart-column"
      ))
    }
    h <- 560
    if (identical(input$opt_box_view_mode, "individual")) {
      h <- 520
    } else {
      df <- get_viz_benchmark_data()
      if (!identical(input$sel_box_cat_group, "all") && !is.null(df)) {
        sub_df <- df[df$Category == input$sel_box_cat_group, , drop = FALSE]
        n_m <- length(unique(sub_df$Metric))
        h <- max(450, min(1400, n_m * 42 + 160))
      } else {
        h <- 780
      }
    }
    plotOutput("plot_metric_box", height = paste0(h, "px"))
  })
  output$plot_metric_box <- renderPlot({ metric_box_reactive() })
  output$download_box_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_boxplots_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, metric_box_reactive(), width = 11, height = 7, dpi = 600) }
  )
  output$download_box_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_boxplots_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 11, height = 7); print(metric_box_reactive()); grDevices::dev.off() }
  )
  output$download_box_cat_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_category_metrics_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, metric_box_reactive(), width = 12, height = 8.5, dpi = 600) }
  )
  output$download_box_cat_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_category_metrics_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 12, height = 8.5); print(metric_box_reactive()); grDevices::dev.off() }
  )
  
  # 5. Metric Heatmap
  output$ui_plot_metric_heat <- renderUI({
    if (rv_plot_triggers$heat == 0) {
      return(render_plot_placeholder(
        title = "Performance Heatmap Ready to Plot",
        desc = "Select category filter and dimensions, then click below to confirm and render the performance matrix.",
        button_id = "btn_plot_heat_inline",
        btn_label = "Confirm & Plot Heatmap",
        icon_name = "table-cells"
      ))
    }
    w <- if (!is.null(input$sld_heat_width)) input$sld_heat_width else 880
    h <- if (!is.null(input$sld_heat_height)) input$sld_heat_height else 1100
    plotOutput("plot_metric_heat", width = paste0(w, "px"), height = paste0(h, "px"))
  })
  
  metric_heat_reactive <- reactive({
    req(rv$benchmark_df)
    req(rv_plot_triggers$heat > 0)
    df <- get_viz_benchmark_data(dataset_choice = input$sel_heat_dataset)
    withProgress(message = "Rendering Performance Heatmap", detail = "Formatting metric matrix...", {
      cat_sel <- if (identical(input$sel_heat_cat, "all")) NULL else input$sel_heat_cat
      b_size <- if (is.null(cat_sel)) 9.5 else 11
      scale_type <- if (!is.null(input$sel_heat_score_type) && identical(input$sel_heat_score_type, "raw")) "raw" else "relative"
      plot_metric_heatmap(
        benchmark_data = df,
        category = cat_sel,
        scale_fill = scale_type,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        base_size = b_size
      )
    })
  })
  output$plot_metric_heat <- renderPlot({ metric_heat_reactive() })
  output$download_heat_jpeg <- downloadHandler(
    filename = function() {
      cat_suffix <- if (identical(input$sel_heat_cat, "all")) "all_categories" else gsub("[^A-Za-z0-9_-]", "_", input$sel_heat_cat)
      paste0("scSimEval_metric_heatmap_", cat_suffix, "_", Sys.Date(), ".jpeg")
    },
    content = function(file) {
      is_all <- is.null(input$sel_heat_cat) || identical(input$sel_heat_cat, "all")
      w <- if (is_all) 9.5 else 9.0
      h <- if (is_all) 15.0 else 7.5
      export_single_jpeg(file, metric_heat_reactive(), width = w, height = h, dpi = 600)
    }
  )
  output$download_heat_pdf <- downloadHandler(
    filename = function() {
      cat_suffix <- if (identical(input$sel_heat_cat, "all")) "all_categories" else gsub("[^A-Za-z0-9_-]", "_", input$sel_heat_cat)
      paste0("scSimEval_metric_heatmap_", cat_suffix, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      is_all <- is.null(input$sel_heat_cat) || identical(input$sel_heat_cat, "all")
      w <- if (is_all) 9.5 else 9.0
      h <- if (is_all) 15.0 else 7.5
      grDevices::pdf(file, width = w, height = h)
      print(metric_heat_reactive())
      grDevices::dev.off()
    }
  )
  
  # 6. PCA Ordination
  observeEvent(input$sel_pca_panel, {
    if (identical(input$sel_pca_panel, "both") || identical(input$sel_pca_panel, "all")) {
      updateSliderInput(session, "sld_pca_height", value = 1100)
    } else {
      updateSliderInput(session, "sld_pca_height", value = 600)
    }
  }, ignoreInit = TRUE)
  
  output$ui_plot_metric_pca <- renderUI({
    if (rv_plot_triggers$pca == 0) {
      return(render_plot_placeholder(
        title = "PCA Ordination Ready to Plot",
        desc = "Configure category filters and PCA view panels, then click below to confirm and compute principal component ordinations.",
        button_id = "btn_plot_pca_inline",
        btn_label = "Confirm & Plot PCA",
        icon_name = "compass"
      ))
    }
    w <- if (!is.null(input$sld_pca_width)) input$sld_pca_width else 950
    h <- if (!is.null(input$sld_pca_height)) input$sld_pca_height else 1100
    plotOutput("plot_metric_pca_out", width = paste0(w, "px"), height = paste0(h, "px"))
  })
  
  metric_pca_reactive <- reactive({
    req(rv$benchmark_df)
    req(rv_plot_triggers$pca > 0)
    df <- get_viz_benchmark_data(dataset_choice = input$sel_pca_dataset)
    n_sims <- length(unique(df$Method))
    if (n_sims < 3) {
      validate(need(FALSE, sprintf("PCA ordination requires at least 3 distinct simulators (currently evaluated: %d). Please evaluate additional methods or select a multi-simulator dataset.", n_sims)))
    }
    withProgress(message = "Computing PCA Ordination", detail = "Decomposing variance and loadings...", {
      cat_sel <- if (identical(input$sel_pca_cat, "all")) NULL else input$sel_pca_cat
      n_top <- if (!is.null(input$num_pca_top_metrics) && is.finite(input$num_pca_top_metrics)) input$num_pca_top_metrics else 14
      p_panel <- if (!is.null(input$sel_pca_panel) && input$sel_pca_panel %in% c("both", "methods", "loadings")) {
        input$sel_pca_panel
      } else if (identical(input$sel_pca_panel, "methods") || identical(input$sel_pca_panel, "scores") || identical(input$sel_pca_panel, "simulators")) {
        "methods"
      } else if (identical(input$sel_pca_panel, "loadings") || identical(input$sel_pca_panel, "scree")) {
        "loadings"
      } else {
        "both"
      }
      tryCatch({
        plot_metric_pca(
          benchmark_data = df,
          category = cat_sel,
          panel = p_panel,
          top_n_loadings = n_top,
          base_size = 12
        )
      }, error = function(e) {
        validate(need(FALSE, paste0("Unable to render PCA plot: ", conditionMessage(e))))
      })
    })
  })
  output$plot_metric_pca_out <- renderPlot({ metric_pca_reactive() })
  output$download_pca_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_pca_", Sys.Date(), ".jpeg") },
    content = function(file) {
      w <- if (!is.null(input$sld_pca_width)) input$sld_pca_width / 90 else 10.5
      h <- if (!is.null(input$sld_pca_height)) input$sld_pca_height / 90 else (if (identical(input$sel_pca_panel, "both")) 12.5 else 7.5)
      export_single_jpeg(file, metric_pca_reactive(), width = w, height = h, dpi = 600)
    }
  )
  output$download_pca_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_pca_", Sys.Date(), ".pdf") },
    content = function(file) {
      w <- if (!is.null(input$sld_pca_width)) input$sld_pca_width / 90 else 10.5
      h <- if (!is.null(input$sld_pca_height)) input$sld_pca_height / 90 else (if (identical(input$sel_pca_panel, "both")) 12.5 else 7.5)
      grDevices::pdf(file, width = w, height = h)
      print(metric_pca_reactive())
      grDevices::dev.off()
    }
  )
  
  # 7. MDS Metric Space
  output$ui_plot_metric_mds <- renderUI({
    if (rv_plot_triggers$mds == 0) {
      return(render_plot_placeholder(
        title = "MDS Distance Projection Ready",
        desc = "Configure distance metric comparison options, then click below to confirm and compute multi-dimensional scaling projection.",
        button_id = "btn_plot_mds_inline",
        btn_label = "Confirm & Plot MDS",
        icon_name = "circle-nodes"
      ))
    }
    w <- if (!is.null(input$sld_mds_width)) input$sld_mds_width else 900
    h <- if (!is.null(input$sld_mds_height)) input$sld_mds_height else 620
    plotOutput("plot_metric_mds_out", width = paste0(w, "px"), height = paste0(h, "px"))
  })
  
  metric_mds_reactive <- reactive({
    req(rv$benchmark_df)
    req(rv_plot_triggers$mds > 0)
    df <- get_viz_benchmark_data(dataset_choice = input$sel_mds_dataset)
    n_sims <- length(unique(df$Method))
    if (n_sims < 3) {
      validate(need(FALSE, sprintf("MDS distance ordination requires at least 3 distinct simulators (currently evaluated: %d). Please evaluate additional methods or select a multi-simulator dataset.", n_sims)))
    }
    withProgress(message = "Computing MDS Metric Space", detail = "Calculating Euclidean distance coordinates...", {
      cat_sel <- if (identical(input$sel_mds_cat, "all")) NULL else input$sel_mds_cat
      ord_by <- if (!is.null(input$sel_mds_by) && input$sel_mds_by %in% c("summaries", "metrics")) {
        "summaries"
      } else {
        "simulators"
      }
      plot_metric_mds(
        benchmark_data = df,
        category = cat_sel,
        ordination_by = ord_by,
        base_size = 13
      )
    })
  })
  output$plot_metric_mds_out <- renderPlot({ metric_mds_reactive() })
  output$download_mds_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_mds_", Sys.Date(), ".jpeg") },
    content = function(file) {
      w <- if (!is.null(input$sld_mds_width)) input$sld_mds_width / 80 else 12
      h <- if (!is.null(input$sld_mds_height)) input$sld_mds_height / 80 else 7.5
      export_single_jpeg(file, metric_mds_reactive(), width = w, height = h, dpi = 600)
    }
  )
  output$download_mds_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_metric_mds_", Sys.Date(), ".pdf") },
    content = function(file) {
      w <- if (!is.null(input$sld_mds_width)) input$sld_mds_width / 80 else 12
      h <- if (!is.null(input$sld_mds_height)) input$sld_mds_height / 80 else 7.5
      grDevices::pdf(file, width = w, height = h)
      print(metric_mds_reactive())
      grDevices::dev.off()
    }
  )
  
  # 8. Cell Embeddings (UMAP, t-SNE, PCA) & Quantitative Quality Metrics
  get_clean_emb_sim_names <- function(cur_dset = 'Dataset 1') {
    if (!is.null(rv$benchmark_df) && "Dataset" %in% colnames(rv$benchmark_df)) {
      sub_bm <- rv$benchmark_df[rv$benchmark_df$Dataset == cur_dset, ]
      if (nrow(sub_bm) > 0) return(unique(sub_bm$Method))
    }
    if (!is.null(rv$methods) && length(rv$methods) > 0) {
      return(rv$methods)
    }
    if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) {
      all_names <- names(rv$sim_matrices)
      clean <- gsub("(_RNA|_ATAC)?(_Dataset.*| \\(Dataset.*\\)|\\[.*\\] )", "", all_names)
      clean <- gsub("^(\\[.*\\]\\s*Simulated:\\s*|Simulated:\\s*)", "", clean)
      clean <- unique(clean[clean != ""])
      if (length(clean) > 0) return(clean)
    }
    return('Simulated')
  }

  # 8. Cell Embeddings (UMAP, t-SNE, PCA) - Collapsible Simulators Selector
  output$ui_emb_sim_picker <- renderUI({
    cur_dset <- input$sel_emb_dataset %||% (input$sel_viz_dataset_view %||% 'Dataset 1')
    if (identical(cur_dset, 'all_avg')) cur_dset <- 'Dataset 1'
    
    sim_names <- get_clean_emb_sim_names(cur_dset)
    sel_methods <- isolate(input$chk_emb_methods)
    if (is.null(sel_methods) || !all(sel_methods %in% sim_names)) {
      sel_methods <- sim_names
    }
    
    tags$details(
      id = 'details_emb_sim_picker',
      class = 'mb-3',
      style = 'background: #F8FAFC; border: 1px solid #CBD5E1; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.04);',
      tags$summary(
        style = 'padding: 8px 12px; background: #EEF2F6; border-radius: 7px; cursor: pointer; list-style: none; font-size: 0.86rem; color: #1B4F72; display: flex; justify-content: space-between; align-items: center;',
        div(
          class = 'd-flex align-items-center gap-2',
          icon('sliders', class = 'text-primary'),
          tags$b('Choose Simulators to Include'),
          span(class = 'badge bg-primary text-white', sprintf('%d / %d', length(sel_methods), length(sim_names)))
        ),
        tags$small(class = 'text-muted', 'Click to hide / expand ▼')
      ),
      div(
        style = 'padding: 10px 12px; background: #FFFFFF; border-top: 1px solid #E2E8F0; border-radius: 0 0 7px 7px;',
        div(
          class = 'd-flex justify-content-between align-items-center mb-2 pb-1 border-bottom',
          tags$span('Simulators:', style = 'font-size: 0.78rem; font-weight: 600; color: #475569;'),
          div(
            actionLink('link_emb_select_all', 'Select All', style = 'font-size: 0.75rem; font-weight: 600; margin-right: 8px;'),
            actionLink('link_emb_clear_all', 'Clear', style = 'font-size: 0.75rem; font-weight: 600; color: #C0392B;')
          )
        ),
        div(
          style = 'max-height: 180px; overflow-y: auto; padding-right: 4px; scrollbar-width: thin;',
          checkboxGroupInput(
            'chk_emb_methods', label = NULL,
            choices = sim_names,
            selected = sel_methods
          )
        )
      )
    )
  })

  observeEvent(input$link_emb_select_all, {
    cur_dset <- input$sel_emb_dataset %||% (input$sel_viz_dataset_view %||% 'Dataset 1')
    if (identical(cur_dset, 'all_avg')) cur_dset <- 'Dataset 1'
    sim_names <- get_clean_emb_sim_names(cur_dset)
    updateCheckboxGroupInput(session, 'chk_emb_methods', selected = sim_names)
  })
  observeEvent(input$link_emb_clear_all, {
    updateCheckboxGroupInput(session, 'chk_emb_methods', selected = character(0))
  })

  output$ui_emb_single_sim <- renderUI({
    req(input$sel_emb_layout == "side_by_side")
    cur_dset <- input$sel_emb_dataset %||% (input$sel_viz_dataset_view %||% 'Dataset 1')
    if (identical(cur_dset, 'all_avg')) cur_dset <- 'Dataset 1'
    sim_names <- get_clean_emb_sim_names(cur_dset)
    if (rv_plot_triggers$emb > 0) {
      emb_data <- tryCatch(reactive_cell_embeddings(), error = function(e) NULL)
      if (!is.null(emb_data) && "Dataset" %in% colnames(emb_data)) {
        emb_sims <- setdiff(unique(as.character(emb_data$Dataset)), "Reference")
        if (length(emb_sims) > 0) sim_names <- emb_sims
      }
    }
    prev_sel <- isolate(input$sel_emb_single_sim)
    sel_val <- if (!is.null(prev_sel) && prev_sel %in% sim_names) prev_sel else sim_names[1]
    selectInput("sel_emb_single_sim", "Select Simulator to Compare:", choices = sim_names, selected = sel_val)
  })

  observe({
    req(input$sel_emb_layout == "side_by_side")
    cur_dset <- input$sel_emb_dataset %||% (input$sel_viz_dataset_view %||% 'Dataset 1')
    if (identical(cur_dset, 'all_avg')) cur_dset <- 'Dataset 1'
    clean_sims <- get_clean_emb_sim_names(cur_dset)
    if (rv_plot_triggers$emb > 0) {
      emb_data <- tryCatch(reactive_cell_embeddings(), error = function(e) NULL)
      if (!is.null(emb_data) && "Dataset" %in% colnames(emb_data)) {
        emb_sims <- setdiff(unique(as.character(emb_data$Dataset)), "Reference")
        if (length(emb_sims) > 0) clean_sims <- emb_sims
      }
    }
    if (!is.null(clean_sims) && length(clean_sims) > 0) {
      cur_val <- isolate(input$sel_emb_single_sim)
      sel_val <- if (!is.null(cur_val) && cur_val %in% clean_sims) cur_val else clean_sims[1]
      updateSelectInput(session, 'sel_emb_single_sim', choices = clean_sims, selected = sel_val)
    }
  })
  
  # Reactive cell embeddings calculation
  reactive_cell_embeddings <- reactive({
    cur_dset <- input$sel_emb_dataset %||% (input$sel_viz_dataset_view %||% "Dataset 1")
    if (identical(cur_dset, "all_avg")) cur_dset <- "Dataset 1"
    cur_mod <- input$sel_emb_modality %||% "RNA"
    cur_dset_clean <- sub(" \\(.*\\)$", "", cur_dset)
    
    # Try resolving reference matrix
    ref_mat <- NULL
    ref_try_keys <- c(
      paste0(cur_dset, "_", cur_mod),
      paste0(cur_dset_clean, "_", cur_mod),
      paste0(cur_dset, " (sc", cur_mod, "-seq)"),
      paste0(cur_dset_clean, " (sc", cur_mod, "-seq)"),
      cur_dset,
      cur_dset_clean,
      paste0(cur_dset, "_", tolower(cur_mod)),
      paste0(cur_dset_clean, "_", tolower(cur_mod))
    )
    for (rk in ref_try_keys) {
      if (!is.null(rv$ref_matrices) && rk %in% names(rv$ref_matrices)) { ref_mat <- rv$ref_matrices[[rk]]; break }
      if (!is.null(rv$all_matrices) && paste0("[Reference] ", rk) %in% names(rv$all_matrices)) { ref_mat <- rv$all_matrices[[paste0("[Reference] ", rk)]]; break }
    }
    if (is.null(ref_mat) && !is.null(rv$ref_matrices) && length(rv$ref_matrices) > 0) {
      mod_hits <- names(rv$ref_matrices)[grepl(cur_mod, names(rv$ref_matrices), ignore.case = TRUE)]
      if (length(mod_hits) > 0) {
        dset_hits <- mod_hits[grepl(cur_dset_clean, mod_hits, fixed = TRUE)]
        ref_mat <- if (length(dset_hits) > 0) rv$ref_matrices[[dset_hits[1]]] else rv$ref_matrices[[mod_hits[1]]]
      } else if (cur_dset %in% names(rv$ref_matrices)) {
        ref_mat <- rv$ref_matrices[[cur_dset]]
      } else {
        ref_mat <- rv$ref_matrices[[1]]
      }
    }
    if (is.null(ref_mat)) ref_mat <- rv$toy_ref
    req(ref_mat)
    req(rv_plot_triggers$emb > 0)
    withProgress(message = "Computing Cell Embeddings", detail = "Running low-dimensional projection...", {
      sim_list <- if (!is.null(rv$sim_matrices) && length(rv$sim_matrices) > 0) {
        rv$sim_matrices
      } else if (!is.null(rv$toy_sim)) {
        list("Simulated" = rv$toy_sim)
      } else {
        NULL
      }
      req(sim_list)
      
      sel_methods <- input$chk_emb_methods
      clean_all <- get_clean_emb_sim_names(cur_dset)
      if (is.null(sel_methods) || length(sel_methods) == 0) sel_methods <- clean_all
      
      # Map chosen simulators to count matrices for the active dataset & modality
      sim_list_filtered <- list()
      for (m in sel_methods) {
        matched_mat <- NULL
        candidate_keys <- c(
          paste0(m, "_", cur_mod, "_", cur_dset),
          paste0(m, "_", cur_mod, "_", cur_dset_clean),
          paste0(m, "_", tolower(cur_mod), "_", cur_dset),
          paste0(m, "_", tolower(cur_mod), "_", cur_dset_clean),
          paste0(m, "_", cur_dset),
          paste0(m, "_", cur_dset_clean),
          paste0(m, "_", cur_mod),
          paste0(m, " (", cur_dset, ")"),
          paste0(m, " (", cur_dset_clean, ")"),
          paste0("Simulated: ", m, " (", cur_dset, ")"),
          paste0("[Simulated] ", m, " ", cur_mod, " (", cur_dset_clean, ")"),
          paste0("[Simulated] ", m, "_", cur_dset),
          paste0("[", cur_dset, "] Simulated: ", m),
          paste0("[", cur_dset, "] Simulated ", cur_mod, ": ", m),
          m,
          paste0("Simulated: ", m)
        )
        for (ck in candidate_keys) {
          if (!is.null(sim_list) && ck %in% names(sim_list)) { matched_mat <- sim_list[[ck]]; break }
          if (!is.null(rv$all_matrices) && ck %in% names(rv$all_matrices)) { matched_mat <- rv$all_matrices[[ck]]; break }
        }
        if (is.null(matched_mat) && !is.null(sim_list)) {
          m_names <- names(sim_list)
          hit <- m_names[grepl(paste0("^", m, "(_| |\\()"), m_names)]
          if (length(hit) > 0) {
            mod_hit <- hit[grepl(cur_mod, hit, ignore.case = TRUE)]
            matched_mat <- if (length(mod_hit) > 0) sim_list[[mod_hit[1]]] else sim_list[[hit[1]]]
          }
        }
        if (!is.null(matched_mat)) {
          sim_list_filtered[[m]] <- matched_mat
        }
      }
      if (length(sim_list_filtered) == 0) {
        sim_list_filtered <- sim_list
        clean_names <- sub("(_RNA|_ATAC)?(_Dataset.*| \\(Dataset.*\\)|\\[.*\\] )", "", names(sim_list_filtered))
        clean_names <- sub("^(\\[.*\\]\\s*Simulated:\\s*|Simulated:\\s*)", "", clean_names)
        names(sim_list_filtered) <- clean_names
      }
      
      find_matching_vector <- function(map_or_val, dset_c, dset_raw, mod_c, n_exp) {
        if (is.null(map_or_val)) return(NULL)
        if (is.list(map_or_val) && !is.data.frame(map_or_val) && !is.factor(map_or_val)) {
          cand_keys <- c(
            paste0(dset_c, "_", mod_c),
            paste0(dset_raw, "_", mod_c),
            paste0(dset_c, "_", tolower(mod_c)),
            paste0(dset_raw, "_", tolower(mod_c)),
            paste0(dset_c, " (sc", mod_c, "-seq)"),
            paste0(dset_raw, " (sc", mod_c, "-seq)"),
            dset_c,
            dset_raw,
            mod_c,
            tolower(mod_c)
          )
          for (k in cand_keys) {
            if (k %in% names(map_or_val)) {
              res <- find_matching_vector(map_or_val[[k]], dset_c, dset_raw, mod_c, n_exp)
              if (!is.null(res)) return(res)
            }
          }
          for (k in cand_keys) {
            m_idx <- match(tolower(k), tolower(names(map_or_val)))
            if (!is.na(m_idx)) {
              res <- find_matching_vector(map_or_val[[m_idx]], dset_c, dset_raw, mod_c, n_exp)
              if (!is.null(res)) return(res)
            }
          }
          for (item in map_or_val) {
            res <- find_matching_vector(item, dset_c, dset_raw, mod_c, n_exp)
            if (!is.null(res)) return(res)
          }
        }
        if (is.data.frame(map_or_val)) {
          v <- if (!is.null(map_or_val$CellType)) map_or_val$CellType else if (!is.null(map_or_val$Cell_Type)) map_or_val$Cell_Type else if (!is.null(map_or_val$Batch)) map_or_val$Batch else map_or_val[[1]]
          if (length(v) == n_exp) return(as.character(v))
        }
        if (length(map_or_val) == n_exp) return(as.character(map_or_val))
        return(NULL)
      }
      
      ct_vec <- find_matching_vector(rv$cell_types_map, cur_dset_clean, cur_dset, cur_mod, ncol(ref_mat))
      if (is.null(ct_vec)) {
        ct_vec <- find_matching_vector(rv$cell_types, cur_dset_clean, cur_dset, cur_mod, ncol(ref_mat))
      }
      bt_vec <- find_matching_vector(rv$batches_map, cur_dset_clean, cur_dset, cur_mod, ncol(ref_mat))
      if (is.null(bt_vec)) {
        bt_vec <- find_matching_vector(rv$batch, cur_dset_clean, cur_dset, cur_mod, ncol(ref_mat))
      }
      
      n_pcs <- if (!is.null(input$sld_emb_pcs)) input$sld_emb_pcs else 20
      perp <- if (!is.null(input$sld_emb_perp)) input$sld_emb_perp else 15
      n_nei <- if (!is.null(input$sld_emb_neighbors)) input$sld_emb_neighbors else 15
      
      compute_dataset_embeddings(
        reference = ref_mat,
        simulated = sim_list_filtered,
        reduction = input$sel_emb_reduction,
        n_pcs = n_pcs,
        perplexity = perp,
        n_neighbors = n_nei,
        cell_types = ct_vec,
        batch = bt_vec
      )
    })
  })
  
  output$ui_plot_cell_embeddings <- renderUI({
    if (rv_plot_triggers$emb == 0) {
      return(render_plot_placeholder(
        title = "Cell Embeddings Ready to Compute",
        desc = "Select projection method (UMAP, t-SNE, PCA) and parameters, then click below to calculate dimensionality reductions.",
        button_id = "btn_plot_emb_inline",
        btn_label = "Compute & Render Embeddings",
        icon_name = "project-diagram"
      ))
    }
    # Embeddings triggered on permission click
    plotOutput("plot_cell_embeddings", height = "650px")
  })
  
  cell_embeddings_plot_reactive <- reactive({
    req(rv_plot_triggers$emb > 0)
    emb_data <- reactive_cell_embeddings()
    req(emb_data)
    lay <- if (identical(input$sel_emb_layout, "side_by_side")) "side_by_side" else "facet"
    pt_sz <- if (!is.null(input$sld_emb_pt_size)) input$sld_emb_pt_size else 1.0
    alpha_val <- if (!is.null(input$sld_emb_alpha)) input$sld_emb_alpha else 0.8
    
    single_sim <- input$sel_emb_single_sim
    if (is.null(single_sim) || !nzchar(trimws(single_sim))) {
      all_sims <- setdiff(unique(as.character(emb_data$Dataset)), "Reference")
      single_sim <- if (length(all_sims) > 0) all_sims[1] else NULL
    }
    
    plot_dataset_embeddings(
      embedding_data = emb_data,
      reduction = input$sel_emb_reduction,
      layout = lay,
      color_by = input$sel_emb_color,
      selected_methods = if (identical(lay, "side_by_side")) single_sim else input$chk_emb_methods,
      pt_size = pt_sz,
      alpha = alpha_val,
      base_size = 12
    )
  })
  output$plot_cell_embeddings <- renderPlot({ cell_embeddings_plot_reactive() })
  output$download_emb_jpeg <- downloadHandler(
    filename = function() { paste0("scSimEval_cell_embeddings_", input$sel_emb_reduction, "_", Sys.Date(), ".jpeg") },
    content = function(file) { export_single_jpeg(file, cell_embeddings_plot_reactive(), width = 13, height = 8.5, dpi = 600) }
  )
  output$download_emb_pdf <- downloadHandler(
    filename = function() { paste0("scSimEval_cell_embeddings_", input$sel_emb_reduction, "_", Sys.Date(), ".pdf") },
    content = function(file) { grDevices::pdf(file, width = 13, height = 8.5); print(cell_embeddings_plot_reactive()); grDevices::dev.off() }
  )
  output$table_emb_quality_metrics <- renderDT({
    req(rv_plot_triggers$emb > 0)
    emb_data <- reactive_cell_embeddings()
    req(emb_data)
    df_metrics <- compute_embedding_quality_metrics(emb_data)
    datatable(df_metrics, options = list(dom = "t", pageLength = 10), rownames = FALSE, class = "compact stripe")
  })
  # Tab 5: Download Results (Data Tables, Excel, CSV, RDS, and Consolidated ZIP)
  # ----------------------------------------------------------------------------
  
  output$download_excel <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".xlsx") },
    content = function(file) {
      req(rv$benchmark_df)
      export_excel_workbook(file, rv$benchmark_df, leaderboard_reactive(), rv$dataset_summary_df)
    }
  )
  
  output$download_csv <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".csv") },
    content = function(file) {
      req(rv$benchmark_df)
      utils::write.csv(rv$benchmark_df, file, row.names = FALSE)
    }
  )
  
  output$download_txt <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".txt") },
    content = function(file) {
      req(rv$benchmark_df)
      utils::write.table(rv$benchmark_df, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  output$download_rds <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_results_", Sys.Date(), ".rds") },
    content = function(file) {
      req(rv$benchmark_df)
      saveRDS(list(
        benchmark_summary_table = rv$benchmark_df,
        methods = rv$methods,
        method_rankings = leaderboard_reactive(),
        dataset_summary = rv$dataset_summary_df
      ), file)
    }
  )
  
  output$download_complete_zip <- downloadHandler(
    filename = function() { paste0("scSimEval_benchmark_data_", Sys.Date(), ".zip") },
    content = function(file) {
      req(rv$benchmark_df)
      
      tmp_dir <- file.path(tempdir(), paste0("scSimEval_data_", as.integer(Sys.time())))
      dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
      
      leader_df <- leaderboard_reactive()
      
      # 1. Excel Workbook (.xlsx)
      export_excel_workbook(file.path(tmp_dir, "scSimEval_benchmark_results.xlsx"), rv$benchmark_df, leader_df, rv$dataset_summary_df)
      
      # 2. Master CSV Table (.csv)
      utils::write.csv(rv$benchmark_df, file.path(tmp_dir, "scSimEval_benchmark_results.csv"), row.names = FALSE)
      
      # 3. Master TSV Table (.tsv)
      utils::write.table(rv$benchmark_df, file.path(tmp_dir, "scSimEval_benchmark_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
      
      # 4. Method Rankings Leaderboard (.csv)
      if (!is.null(leader_df)) {
        utils::write.csv(leader_df, file.path(tmp_dir, "scSimEval_method_rankings.csv"), row.names = FALSE)
      }
      
      # 5. Dataset Properties Summary (.csv)
      if (!is.null(rv$dataset_summary_df)) {
        utils::write.csv(rv$dataset_summary_df, file.path(tmp_dir, "scSimEval_dataset_properties.csv"), row.names = FALSE)
      }
      
      # 6. Serialized R Data Object (.rds)
      saveRDS(list(
        benchmark_summary_table = rv$benchmark_df,
        methods = rv$methods,
        method_rankings = leader_df,
        dataset_summary = rv$dataset_summary_df
      ), file.path(tmp_dir, "scSimEval_benchmark_data.rds"))
      
      # 7. Informative README
      readme_lines <- c(
        "scSimEval Benchmark Results Data Archive",
        "========================================",
        paste0("Generated on: ", Sys.time()),
        paste0("Source: ", rv$source_name %||% "scSimEval Shiny Evaluation"),
        "",
        "Files included in this archive:",
        "  - scSimEval_benchmark_results.xlsx : Full results workbook with rankings & category summaries",
        "  - scSimEval_benchmark_results.csv  : Master benchmark scores table in CSV format",
        "  - scSimEval_benchmark_results.tsv  : Master benchmark scores table in TSV format",
        if (!is.null(leader_df)) "  - scSimEval_method_rankings.csv    : Overall simulator leaderboard rankings" else NULL,
        if (!is.null(rv$dataset_summary_df)) "  - scSimEval_dataset_properties.csv : Dataset metadata and quality metrics" else NULL,
        "  - scSimEval_benchmark_data.rds     : Complete R list object for downstream scripting in R",
        "",
        "Note: Individual high-resolution publication figures (JPEG 600 DPI, PDF) can be customized",
        "and downloaded directly from the Comparative Bubble Matrix and Visualizations tabs."
      )
      writeLines(stats::na.omit(readme_lines), file.path(tmp_dir, "README_scSimEval.txt"))
      
      zip_files <- list.files(tmp_dir, full.names = FALSE, recursive = TRUE)
      zip::zip(file, files = zip_files, root = tmp_dir)
      unlink(tmp_dir, recursive = TRUE)
    }
  )
  
  # Unified format downloader based on sel_download_format
  # Dedicated Master Table Download Handlers inside the Table Box
  output$download_master_table_excel <- downloadHandler(
    filename = function() {
      paste0("scSimEval_master_benchmark_scores_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(rv$benchmark_df)
      df_export <- rv$benchmark_df
      if (!is.null(input$table_master_export_rows_all) &&
          length(input$table_master_export_rows_all) > 0 &&
          length(input$table_master_export_rows_all) < nrow(df_export)) {
        df_export <- df_export[input$table_master_export_rows_all, , drop = FALSE]
      }
      export_excel_workbook(file, df_export, leaderboard_reactive(), rv$dataset_summary_df)
    }
  )

  output$download_master_table_csv <- downloadHandler(
    filename = function() {
      paste0("scSimEval_master_benchmark_scores_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$benchmark_df)
      df_export <- rv$benchmark_df
      if (!is.null(input$table_master_export_rows_all) &&
          length(input$table_master_export_rows_all) > 0 &&
          length(input$table_master_export_rows_all) < nrow(df_export)) {
        df_export <- df_export[input$table_master_export_rows_all, , drop = FALSE]
      }
      utils::write.csv(df_export, file, row.names = FALSE)
    }
  )

  # Benchmark Results Session Export (Exclusively .rds for reproduction)
  output$download_benchmark_data <- downloadHandler(
    filename = function() {
      paste0("scSimEval_benchmark_session_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(rv$benchmark_df)
      fmt <- input$sel_download_format %||% "xlsx"
      
      # Export active filtered rows if user has filtered the table in DT
      df_export <- rv$benchmark_df
      if (!is.null(input$table_master_export_rows_all) && 
          length(input$table_master_export_rows_all) > 0 && 
          length(input$table_master_export_rows_all) < nrow(df_export)) {
        df_export <- df_export[input$table_master_export_rows_all, , drop = FALSE]
      }
      
      # Retain both raw evaluation metric values and normalized scores in exported table
      raw_val <- if ("Score_Raw" %in% colnames(df_export)) {
        df_export$Score_Raw
      } else if ("Value" %in% colnames(df_export)) {
        df_export$Value
      } else if ("Raw_Value" %in% colnames(df_export)) {
        df_export$Raw_Value
      } else {
        df_export$Score
      }
      norm_val <- if ("Score" %in% colnames(df_export)) {
        df_export$Score
      } else if ("Score_Norm" %in% colnames(df_export)) {
        df_export$Score_Norm
      } else {
        raw_val
      }
      df_export$Score_Norm <- NULL
      df_export$Value <- NULL
      df_export$Raw_Value <- NULL
      df_export$Score_Raw <- round(as.numeric(raw_val), 4)
      df_export$Score <- round(as.numeric(norm_val), 4)
      pref_cols <- c("Dataset", "Method", "Category", intersect("Property", colnames(df_export)), "Metric", "Score_Raw", "Score")
      col_order <- unique(c(pref_cols, colnames(df_export)))
      col_order <- intersect(col_order, colnames(df_export))
      df_export <- df_export[, col_order, drop = FALSE]
      
      if (fmt == "xlsx") {
        export_excel_workbook(file, df_export, leaderboard_reactive(), rv$dataset_summary_df)
      } else if (fmt == "csv") {
        utils::write.csv(df_export, file, row.names = FALSE)
      } else if (fmt == "tsv") {
        utils::write.table(df_export, file, sep = "\t", row.names = FALSE, quote = FALSE)
      } else if (fmt == "rds") {
        saveRDS(list(
          benchmark_summary_table = df_export,
          methods = rv$methods,
          method_rankings = leaderboard_reactive(),
          dataset_summary = rv$dataset_summary_df
        ), file)
      } else if (fmt == "zip") {
        tmp_dir <- file.path(tempdir(), paste0("scSimEval_data_", as.integer(Sys.time())))
        dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
        leader_df <- leaderboard_reactive()
        
        export_excel_workbook(file.path(tmp_dir, "scSimEval_benchmark_results.xlsx"), df_export, leader_df, rv$dataset_summary_df)
        utils::write.csv(df_export, file.path(tmp_dir, "scSimEval_benchmark_results.csv"), row.names = FALSE)
        utils::write.table(df_export, file.path(tmp_dir, "scSimEval_benchmark_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
        if (!is.null(leader_df)) utils::write.csv(leader_df, file.path(tmp_dir, "scSimEval_method_rankings.csv"), row.names = FALSE)
        if (!is.null(rv$dataset_summary_df)) utils::write.csv(rv$dataset_summary_df, file.path(tmp_dir, "scSimEval_dataset_properties.csv"), row.names = FALSE)
        saveRDS(list(
          benchmark_summary_table = df_export,
          methods = rv$methods,
          method_rankings = leader_df,
          dataset_summary = rv$dataset_summary_df
        ), file.path(tmp_dir, "scSimEval_benchmark_data.rds"))
        
        readme_lines <- c(
          "scSimEval Benchmark Results Data Archive",
          "========================================",
          paste0("Generated on: ", Sys.time()),
          paste0("Source: ", rv$source_name %||% "scSimEval Shiny Evaluation"),
          "",
          "Files included in this archive:",
          "  - scSimEval_benchmark_results.xlsx : Full results workbook with rankings & category summaries",
          "  - scSimEval_benchmark_results.csv  : Master benchmark scores table in CSV format",
          "  - scSimEval_benchmark_results.tsv  : Master benchmark scores table in TSV format",
          if (!is.null(leader_df)) "  - scSimEval_method_rankings.csv    : Overall simulator leaderboard rankings" else NULL,
          if (!is.null(rv$dataset_summary_df)) "  - scSimEval_dataset_properties.csv : Dataset metadata and quality metrics" else NULL,
          "  - scSimEval_benchmark_data.rds     : Complete R list object for downstream scripting in R",
          "",
          "Note: Individual high-resolution publication figures (JPEG 600 DPI, PDF) can be customized",
          "and downloaded directly from the Comparative Bubble Matrix and Visualizations tabs."
        )
        writeLines(stats::na.omit(readme_lines), file.path(tmp_dir, "README_scSimEval.txt"))
        
        zip_files <- list.files(tmp_dir, full.names = FALSE, recursive = TRUE)
        zip::zip(file, files = zip_files, root = tmp_dir)
        unlink(tmp_dir, recursive = TRUE)
      }
    }
  )
  
  output$table_master_export <- renderDT({
    req(rv$benchmark_df)
    df <- rv$benchmark_df
    
    # Ensure Dataset column exists and is positioned first for instant dropdown filtering
    if (!"Dataset" %in% colnames(df)) df$Dataset <- "Dataset 1"
    
    # Ensure raw values and normalized scores are both present
    if (!"Score_Raw" %in% colnames(df)) {
      if ("Value" %in% colnames(df)) df$Score_Raw <- df$Value
      else if ("Raw_Value" %in% colnames(df)) df$Score_Raw <- df$Raw_Value
      else df$Score_Raw <- df$Score
    }
    if (!"Score" %in% colnames(df) && "Score_Norm" %in% colnames(df)) {
      df$Score <- df$Score_Norm
    } else if (!"Score" %in% colnames(df)) {
      df$Score <- df$Score_Raw
    }
    
    # Ensure Score is rounded normalized score in [0, 1]
    df$Score <- round(as.numeric(df$Score), 4)
    df$Score_Raw <- round(as.numeric(df$Score_Raw), 4)
    
    # Ensure Property column is present and never empty
    if (!"Property" %in% colnames(df)) {
      df$Property <- df$Metric
    } else {
      df$Property <- ifelse(is.na(df$Property) | df$Property == "", df$Metric, df$Property)
    }
    
    # Remove redundant intermediate columns
    df$Score_Norm <- NULL
    df$Value <- NULL
    df$Raw_Value <- NULL
    
    # Establish clean, logical column ordering:
    # Dataset -> Method -> Category -> Property -> Metric -> Score_Raw -> Score -> Modality
    pref_cols <- c("Dataset", "Method", "Category", "Property", "Metric", "Score_Raw", "Score", intersect("Modality", colnames(df)))
    col_order <- unique(c(pref_cols, colnames(df)))
    col_order <- intersect(col_order, colnames(df))
    df <- df[, col_order, drop = FALSE]
    
    # Convert text columns to factors for automatic DT dropdown filtering
    cat_cols <- c("Dataset", "Method", "Category", "Metric", "Property")
    for (col in cat_cols) {
      if (col %in% colnames(df)) df[[col]] <- as.factor(df[[col]])
    }
    
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
  
  # ============================================================================
  # TAB 8: LOGS SERVER LOGIC & CONSOLE RENDERING
  # ============================================================================
  
  # Reactive expression returning filtered log entries
  filtered_log_entries <- reactive({
    input$btn_refresh_logs
    sc_log_version()
    
    entries <- log_storage$entries
    if (length(entries) == 0) return(list())
    
    # Severity filter
    lvl <- input$log_filter_level %||% "ALL"
    if (lvl != "ALL") {
      entries <- Filter(function(e) isTRUE(e$level == lvl), entries)
    }
    
    # Text search
    q <- trimws(input$log_search_query %||% "")
    if (nzchar(q)) {
      ql <- tolower(q)
      entries <- Filter(function(e) {
        grepl(ql, tolower(e$message), fixed = TRUE) ||
          grepl(ql, tolower(e$component), fixed = TRUE) ||
          (!is.null(e$details) && grepl(ql, tolower(e$details), fixed = TRUE))
      }, entries)
    }
    
    return(entries)
  })
  
  # KPI Metric Output Badges
  output$log_stat_total <- renderUI({
    sc_log_version()
    length(log_storage$entries)
  })
  
  output$log_stat_error <- renderUI({
    sc_log_version()
    sum(vapply(log_storage$entries, function(e) identical(e$level, "ERROR"), logical(1)))
  })
  
  output$log_stat_warning <- renderUI({
    sc_log_version()
    sum(vapply(log_storage$entries, function(e) identical(e$level, "WARNING"), logical(1)))
  })
  
  output$log_stat_note_info <- renderUI({
    sc_log_version()
    sum(vapply(log_storage$entries, function(e) e$level %in% c("NOTE", "INFO", "SUCCESS"), logical(1)))
  })
  
  output$log_console_status_line <- renderUI({
    sc_log_version()
    f_count <- length(filtered_log_entries())
    t_count <- length(log_storage$entries)
    sprintf("Showing %d of %d events", f_count, t_count)
  })
  
  # Render Live Console HTML
  output$ui_log_console <- renderUI({
    entries <- filtered_log_entries()
    
    if (length(entries) == 0) {
      return(
        div(
          style = "text-align: center; padding: 60px 20px; color: #64748B;",
          icon("terminal", style = "font-size: 2.4rem; opacity: 0.4; margin-bottom: 12px;"),
          h6("No Log Events Match Filter Criteria", style = "color: #94A3B8; font-weight: 700;"),
          p("Execute operations in the Data Hub or Visualizations tabs to generate real-time execution logs.", style = "font-size: 0.82rem; margin-bottom: 0;")
        )
      )
    }
    
    rendered_rows <- lapply(entries, function(e) {
      lvl_cls <- paste0("log-lvl-badge log-lvl-", e$level)
      
      div(
        class = "log-entry-row",
        span(class = "log-time", e$time),
        span(class = lvl_cls, e$level),
        span(class = "log-comp", paste0("[", e$component, "]")),
        div(
          class = "log-text",
          span(e$message),
          if (!is.null(e$details) && nzchar(e$details)) {
            div(class = "log-details-block", tags$pre(style = "margin: 0; font-family: inherit; color: inherit; white-space: pre-wrap;", e$details))
          }
        )
      )
    })
    
    autoscroll_script <- if (isTRUE(input$log_autoscroll)) {
      tags$script(HTML("setTimeout(function(){ var el = document.getElementById('sc_log_body'); if(el) { el.scrollTop = el.scrollHeight; } }, 50);"))
    } else {
      NULL
    }
    
    tagList(
      rendered_rows,
      autoscroll_script
    )
  })
  
  # Clear Logs Action
  observeEvent(input$btn_clear_logs, {
    log_storage$entries <- list()
    log_storage$counter <- 0L
    sc_log_version(isolate(sc_log_version()) + 1L)
    sc_log("INFO", "Diagnostic log history cleared by user session.", "Console")
    showNotification("Console logs cleared.", type = "message", duration = 3)
  })
  
  # Refresh Logs Action
  observeEvent(input$btn_refresh_logs, {
    sc_log_version(isolate(sc_log_version()) + 1L)
  })
  
  # Download Logs Handler (.txt)
  output$btn_download_logs <- downloadHandler(
    filename = function() {
      paste0("scSimEval_logs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
    },
    content = function(file) {
      entries <- log_storage$entries
      hdr <- c(
        "================================================================================",
        " scSimEval: Single-Cell Multiomics Simulation Benchmarking Studio",
        " Interactive Execution & Diagnostic Runtime Log Export",
        paste0(" Generated At:    ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste0(" Package Version: scSimEval v0.99.9"),
        paste0(" R Version:       ", R.version.string),
        paste0(" Platform:        ", R.version$platform),
        paste0(" Total Entries:   ", length(entries)),
        "================================================================================",
        ""
      )
      
      body_lines <- if (length(entries) > 0) {
        unlist(lapply(entries, function(e) {
          line <- sprintf("[%s] [%-7s] [%-12s] %s", e$time, e$level, e$component, e$message)
          if (!is.null(e$details) && nzchar(e$details)) {
            c(line, paste0("    Trace/Details: ", e$details))
          } else {
            line
          }
        }))
      } else {
        "[No log entries recorded]"
      }
      
      writeLines(c(hdr, body_lines), con = file)
    }
  )
}

# ==============================================================================
# Run Shiny App
# ==============================================================================
shinyApp(ui = ui, server = server)
