test_that("one-condition sites are retained and classified as not contrastable", {
  skip_if_not_installed("glmmTMB")
  dat <- data.frame(
    site_id=rep(c("s1","s2"), each=3),
    ed=c(2,3,4, 1,2,3), total=rep(10,6),
    conditions=c("A","A","A", "A","B","B")
  )
  fit <- fit_bb_repeated(
    dat, site="site_id", family_order="binomial",
    min_samples=3, quiet=TRUE
  )
  expect_equal(length(fit$fits), 2L)
  expect_equal(fit$fits[["s1"]]$status, "not_contrastable")
  expect_equal(fit$fits[["s1"]]$failure_reason, "fewer than two condition levels")
  expect_equal(fit$fits[["s2"]]$status, "success")
  expect_equal(fit$n_sites_retained, 2L)
  expect_equal(fit$n_sites_attempted, 1L)
  expect_equal(fit$n_sites_successful, 1L)
  expect_equal(fit$n_sites_failed, 0L)
  expect_equal(fit$n_sites_not_contrastable, 1L)
})

test_that("custom contrasts unavailable at one site do not abort other sites", {
  skip_if_not_installed("glmmTMB")
  dat <- data.frame(
    site_id=rep(c("s1","s2"), each=4),
    ed=c(2,3,4,5, 2,3,4,5), total=rep(10,8),
    conditions=c("A","A","B","B", "A","A","C","C")
  )
  fit <- fit_bb_repeated(
    dat, site="site_id", family_order="binomial",
    contrast_pairs=data.frame(A="A", B="B"), quiet=TRUE
  )
  expect_equal(length(fit$fits), 2L)
  expect_equal(fit$fits[["s1"]]$status, "success")
  expect_equal(fit$fits[["s2"]]$status, "not_contrastable")
  expect_true(is.null(fit$fits[["s2"]]$fit))
  expect_false(isTRUE(fit$fits[["s2"]]$fit_attempted))
  expect_match(fit$fits[["s2"]]$failure_reason, "not present")
  expect_equal(fit$n_sites_successful, 1L)
  expect_equal(fit$n_sites_not_contrastable, 1L)
})

test_that("compare_bb_methods distinguishes raw and adjusted significance", {
  x <- tibble::tibble(
    coefficient_diff=c(1, 1), delta_diff=c(0.4, 0.4),
    coefficient_p=c(0.01, 0.04), delta_p=c(0.02, 0.20),
    coefficient_p_adj=c(0.08, 0.08), delta_p_adj=c(0.08, 0.20)
  )
  out <- compare_bb_methods(x)
  expect_true(all(c(
    "coefficient_significant_raw", "delta_significant_raw",
    "significance_agreement_raw", "coefficient_significant_adj",
    "delta_significant_adj", "significance_agreement_adj",
    "significance_agreement"
  ) %in% names(out)))
  expect_equal(out$coefficient_significant_raw, c(TRUE, TRUE))
  expect_equal(out$delta_significant_raw, c(TRUE, FALSE))
  expect_equal(out$significance_agreement_raw, c(TRUE, FALSE))
  expect_equal(out$coefficient_significant_adj, c(FALSE, FALSE))
  expect_equal(out$delta_significant_adj, c(FALSE, FALSE))
  expect_equal(out$significance_agreement_adj, c(TRUE, TRUE))
  expect_equal(out$significance_agreement, out$significance_agreement_adj)
})

test_that("compare_bb_methods falls back to raw significance when adjusted p-values are absent", {
  x <- tibble::tibble(
    coefficient_diff=1, delta_diff=0.4,
    coefficient_p=0.01, delta_p=0.20
  )
  out <- compare_bb_methods(x)
  expect_true("significance_agreement_raw" %in% names(out))
  expect_equal(out$significance_agreement, FALSE)
})
