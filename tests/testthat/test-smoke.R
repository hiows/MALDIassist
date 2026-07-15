test_that("find_peaks_fast detects synthetic peaks", {
  x <- seq(1000, 2000, length.out = 500)
  y <- dnorm(x, mean = 1400, sd = 20) * 100 +
    dnorm(x, mean = 1600, sd = 20) * 80
  spec <- data.frame(mz = x, intensity = y)

  peaks <- find_peaks_fast(spec, hws_peaks = 10)

  expect_s3_class(peaks, "data.frame")
  expect_true(nrow(peaks) >= 2)
  expect_true(all(c("mz", "intensity") %in% names(peaks)))
})

test_that("estimate_significance returns p-values for two groups", {
  set.seed(1)
  n_per_group <- 8
  mat <- rbind(
    matrix(rnorm(n_per_group * 4, mean = 0), nrow = n_per_group),
    matrix(rnorm(n_per_group * 4, mean = 1), nrow = n_per_group)
  )
  colnames(mat) <- paste0("mz_", 1:4)
  group <- rep(c("A", "B"), each = n_per_group)

  out <- estimate_significance(mat, group, adj_method = "BH")

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), ncol(mat))
  expect_true(all(c("feat_names", "pvalue", "adj_pvalue") %in% names(out)))
})

test_that("build_matched_matrix assembles sample-by-marker matrix", {
  peaks_list <- list(
    s1 = data.frame(mz = c(1000, 2000), intensity = c(10, 20)),
    s2 = data.frame(mz = c(1005, 1990), intensity = c(12, 18))
  )

  matched <- build_matched_matrix(
    peaks_list = peaks_list,
    reference_mz = c(1000, 2000),
    hws_match = 20
  )

  expect_type(matched, "list")
  expect_true(is.matrix(matched$detected_matrix))
  expect_equal(nrow(matched$detected_matrix), 2)
  expect_equal(ncol(matched$detected_matrix), 2)
})

test_that("subtract_baseline returns corrected spectrum data", {
  x <- seq(1, 100, length.out = 200)
  y <- sin(x / 5) + 5
  spec <- data.frame(mz = x, intensity = y)

  out <- subtract_baseline(
    data = spec,
    baseline_type = "snip",
    iter_snip = 20
  )

  expect_type(out, "list")
  expect_true(is.data.frame(out$subtracted_data))
  expect_equal(nrow(out$subtracted_data), nrow(spec))
  expect_true(max(out$subtracted_data[, 2]) <= max(spec$intensity))
})
