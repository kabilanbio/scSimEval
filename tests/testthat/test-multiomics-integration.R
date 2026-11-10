test_that("Category VI integration metrics compute accurately on synthetic embeddings", {
  set.seed(42)
  n <- 120
  ct <- sample(c("TypeA", "TypeB", "TypeC"), n, replace = TRUE)
  bt <- sample(c("Batch1", "Batch2"), n, replace = TRUE)
  emb <- matrix(rnorm(n * 10), nrow = n, ncol = 10)
  
  # Inject cell-type separation
  for (i in seq_len(n)) {
    if (ct[i] == "TypeA") emb[i, 1:3] <- emb[i, 1:3] + 4
    if (ct[i] == "TypeB") emb[i, 4:6] <- emb[i, 4:6] + 4
  }
  
  b_asw <- calc_batch_asw(emb, bt, ct)
  expect_true(is.numeric(b_asw) && b_asw >= 0 && b_asw <= 1)
  
  ct_asw <- calc_cell_type_asw(emb, ct)
  expect_true(is.numeric(ct_asw) && ct_asw >= 0 && ct_asw <= 1)
  
  g_conn <- calc_graph_connectivity(emb, ct, k = 15)
  expect_true(is.numeric(g_conn) && g_conn >= 0 && g_conn <= 1)
  
  ari_val <- calc_integration_ari(emb, ct)
  expect_true(is.numeric(ari_val) && ari_val >= 0 && ari_val <= 1)
  
  nmi_val <- calc_integration_nmi(emb, ct)
  expect_true(is.numeric(nmi_val) && nmi_val >= 0 && nmi_val <= 1)
})

test_that("evaluate_multiomics_integration handles multiomics lists correctly", {
  data("example_multiomics", package = "scSimEval", envir = environment())
  m <- example_multiomics
  
  ct <- rep(c("CellType_A", "CellType_B"), each = as.integer(ncol(m$sim_multi[[1]]) / 2))
  bt <- rep(c("Batch_1", "Batch_2"), length.out = ncol(m$sim_multi[[1]]))
  
  res <- evaluate_multiomics_integration(
    ref_multi = m$ref_multi,
    sim_multi = m$sim_multi,
    cell_types = ct,
    batch_info = bt,
    verbose = FALSE
  )
  
  expect_s3_class(res, "scSimEval_integration")
  expect_type(res$metrics, "list")
  expect_equal(length(res$metrics), 5)
  
  expected_metrics <- c("batch_asw", "cell_type_asw", "graph_connectivity", "ari", "nmi")
  expect_true(all(expected_metrics %in% names(res$metrics)))
  
  expect_s3_class(res$integration_summary_table, "data.frame")
  expect_equal(nrow(res$integration_summary_table), 5)
  expect_true(all(res$integration_summary_table$Metric %in% expected_metrics))
  expect_true(all(res$integration_summary_table$Value >= 0 & res$integration_summary_table$Value <= 1))
})

test_that("evaluate_multiomics_accuracy integrates Category VI metrics into master summary", {
  data("example_multiomics", package = "scSimEval", envir = environment())
  m <- example_multiomics
  
  res <- evaluate_multiomics_accuracy(
    ref_multi = m$ref_multi,
    sim_multi = m$sim_multi,
    cell_types = rep(c("A", "B"), each = 50),
    batch_info = rep(c("Rep1", "Rep2"), 50),
    compute_bivariate = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(res$benchmark_summary_table, "data.frame")
  std_tbl <- standardize_benchmark_categories(res$benchmark_summary_table)
  
  cat_vi_rows <- std_tbl[std_tbl$Category == "(VI) Multiomics Integration & Alignment", , drop = FALSE]
  expect_true(nrow(cat_vi_rows) >= 5)
  
  expected_metrics <- c("batch_asw", "cell_type_asw", "graph_connectivity", "ari", "nmi")
  expect_true(all(expected_metrics %in% cat_vi_rows$Metric))
})

