test_that(".normalize_bubble_scores keeps 2D KS distances on a proper min-max scale", {
  df <- data.frame(
    Method = c("Method_A", "Method_B"),
    Metric = rep("Fasano_Franceschini_2D_KS", 2),
    Score_Raw = c(850, 100),
    stringsAsFactors = FALSE
  )

  out <- .normalize_bubble_scores(df)

  expect_true(all(out$Score_Norm >= 0 & out$Score_Norm <= 1))
  expect_equal(out$Score_Norm, c(0, 1), tolerance = 1e-8)
  expect_false(all(out$Score_Norm == 0))
})
