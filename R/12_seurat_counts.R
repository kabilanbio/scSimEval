# Seurat count extraction compatibility helpers
#
# Keep all Seurat count extraction behind these helpers. Seurat v5 Assay5
# stores expression data in layers (including split counts.* layers), not in
# the legacy @counts slot used by Assay.

.seurat_assay_counts <- function(obj, assay_name = NULL) {
  if (!inherits(obj, "Seurat")) {
    stop("`obj` must be a Seurat object to extract assay counts.", call. = FALSE)
  }
  if (is.null(assay_name) || length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
    assay_name <- .seurat_default_assay(obj)
  }
  if (is.null(assay_name) || !assay_name %in% names(obj@assays)) {
    stop("The requested Seurat assay is not present in the object.", call. = FALSE)
  }

  assay_obj <- obj[[assay_name]]
  if (is.null(assay_obj)) {
    stop(sprintf("Seurat assay '%s' could not be retrieved.", assay_name), call. = FALSE)
  }

  # Assay5 may contain one counts layer or multiple split layers (counts.*).
  # Join split layers on a local assay copy; do not modify the caller's object.
  if (inherits(assay_obj, "Assay5")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("SeuratObject is required to read counts from a Seurat v5 Assay5.", call. = FALSE)
    }

    layer_names <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character())
    count_layers <- grep("^counts($|\\.)", layer_names, value = TRUE)
    if (length(count_layers) == 0L) {
      stop(sprintf("Seurat assay '%s' has no raw counts layer.", assay_name), call. = FALSE)
    }

    assay_for_read <- assay_obj
    if (length(count_layers) > 1L || !identical(count_layers, "counts")) {
      assay_for_read <- tryCatch(
        SeuratObject::JoinLayers(assay_obj, layers = count_layers, new = "counts"),
        error = function(e) {
          stop(sprintf("Could not join split counts layers in Seurat assay '%s': %s", assay_name, conditionMessage(e)), call. = FALSE)
        }
      )
    }

    counts <- tryCatch(
      SeuratObject::LayerData(assay_for_read, layer = "counts"),
      error = function(e) {
        tryCatch(
          Seurat::GetAssayData(obj, assay = assay_name, layer = "counts"),
          error = function(e2) NULL
        )
      }
    )
    if (!is.null(counts) && length(dim(counts)) == 2L && nrow(counts) > 0L) {
      return(counts)
    }
    stop(sprintf("Could not retrieve the counts layer from Seurat assay '%s'.", assay_name), call. = FALSE)
  }

  # Seurat v3/v4 Assay objects expose counts through GetAssayData(slot=).
  counts <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay_name, slot = "counts"),
    error = function(e) NULL
  )
  if (!is.null(counts) && length(dim(counts)) == 2L && nrow(counts) > 0L) {
    return(counts)
  }

  # Also support assay implementations that expose the modern layer API.
  counts <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay_name, layer = "counts"),
    error = function(e) NULL
  )
  if (!is.null(counts) && length(dim(counts)) == 2L && nrow(counts) > 0L) {
    return(counts)
  }

  stop(sprintf("Could not retrieve raw counts from Seurat assay '%s'.", assay_name), call. = FALSE)
}

.get_seurat_counts <- function(obj) {
  .seurat_assay_counts(obj, .seurat_default_assay(obj))
}
