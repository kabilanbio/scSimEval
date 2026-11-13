test_that("Seurat v5 count extraction supports Assay5 and split counts layers", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Seurat")

  counts <- matrix(
    seq_len(24), nrow = 6, ncol = 4,
    dimnames = list(paste0("feature", 1:6), paste0("cell", 1:4))
  )
  obj <- SeuratObject::CreateSeuratObject(counts = counts)
  assay <- obj[[SeuratObject::DefaultAssay(obj)]]
  skip_if_not(inherits(assay, "Assay5"), "Requires SeuratObject v5 Assay5")

  expect_equal(scSimEval:::.seurat_assay_counts(obj), counts)
  expect_equal(scSimEval:::.get_seurat_counts(obj), counts)
  expect_equal(scSimEval::validate_sc_object(obj)$counts, counts)

  # Seurat v5 can represent counts as per-sample layers. Extraction should
  # join those layers on a copy and return the original feature-by-cell matrix.
  obj_split <- obj
  obj_split[[SeuratObject::DefaultAssay(obj_split)]] <- split(
    obj_split[[SeuratObject::DefaultAssay(obj_split)]],
    f = factor(c("sample1", "sample1", "sample2", "sample2"))
  )
  expect_equal(scSimEval:::.seurat_assay_counts(obj_split), counts)
  expect_equal(scSimEval::validate_sc_object(obj_split)$counts, counts)
})

test_that("Seurat v5 multiomics count extraction works for both Assay5 modalities", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Seurat")

  rna <- matrix(
    seq_len(24), nrow = 6, ncol = 4,
    dimnames = list(paste0("gene", 1:6), paste0("cell", 1:4))
  )
  atac <- matrix(
    seq_len(20), nrow = 5, ncol = 4,
    dimnames = list(paste0("peak", 1:5), paste0("cell", 1:4))
  )
  obj <- SeuratObject::CreateSeuratObject(counts = rna, assay = "RNA")
  obj[["ATAC"]] <- SeuratObject::CreateAssay5Object(counts = atac)
  skip_if_not(inherits(obj[["RNA"]], "Assay5") && inherits(obj[["ATAC"]], "Assay5"),
              "Requires SeuratObject v5 Assay5")

  validated <- scSimEval::validate_multiomics_input(obj)
  expect_equal(validated$m1, rna)
  expect_equal(validated$m2, atac)
})
