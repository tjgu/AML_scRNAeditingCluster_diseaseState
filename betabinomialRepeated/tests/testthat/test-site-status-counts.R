test_that("unavailable requested contrasts are not attempted and are not successful", {
  skip_if_not_installed("glmmTMB")

  dat <- data.frame(
    site = c(rep("s1", 4), rep("s2", 4)),
    condition = c("A", "B", "A", "B", "A", "C", "A", "C"),
    edited = c(1, 2, 2, 3, 1, 2, 2, 3),
    total = rep(10, 8)
  )

  fit <- fit_bb_repeated(
    data = dat,
    site = "site",
    edited = "edited",
    total = "total",
    condition = "condition",
    contrast_pairs = data.frame(A = "A", B = "B"),
    family_order = c("binomial")
  )

  expect_equal(fit$n_sites_retained, 2L)
  expect_equal(fit$n_sites_attempted, 1L)
  expect_equal(fit$n_sites_successful, 1L)
  expect_equal(fit$n_sites_failed, 0L)
  expect_equal(fit$n_sites_not_contrastable, 1L)

  expect_identical(fit$fits[["s1"]]$status, "success")
  expect_identical(fit$fits[["s2"]]$status, "not_contrastable")
  expect_false(isTRUE(fit$fits[["s2"]]$fit_attempted))
  expect_true(is.null(fit$fits[["s2"]]$fit))
  expect_match(fit$fits[["s2"]]$failure_reason, "not present")
})

test_that("one-condition sites are retained but are not model attempts", {
  skip_if_not_installed("glmmTMB")

  dat <- data.frame(
    site = rep("s1", 4),
    condition = rep("A", 4),
    edited = c(1, 2, 1, 2),
    total = rep(10, 4)
  )

  fit <- fit_bb_repeated(
    data = dat,
    site = "site",
    edited = "edited",
    total = "total",
    condition = "condition",
    family_order = c("binomial")
  )

  expect_equal(fit$n_sites_retained, 1L)
  expect_equal(fit$n_sites_attempted, 0L)
  expect_equal(fit$n_sites_successful, 0L)
  expect_equal(fit$n_sites_failed, 0L)
  expect_equal(fit$n_sites_not_contrastable, 1L)
  expect_identical(fit$fits[["s1"]]$status, "not_contrastable")
  expect_false(isTRUE(fit$fits[["s1"]]$fit_attempted))
})
