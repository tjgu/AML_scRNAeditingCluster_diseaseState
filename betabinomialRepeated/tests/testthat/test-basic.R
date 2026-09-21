test_that("input validation rejects invalid family_order", {
  dat <- data.frame(
    site_id=rep("s1", 4), ed=c(2,3,4,5), total=rep(10,4),
    conditions=c("A","A","B","B")
  )
  expect_error(
    fit_bb_repeated(dat, site="site_id", family_order="invalid"),
    "family_order"
  )
})

test_that("missing categorical identifiers are rejected", {
  dat <- data.frame(
    site_id=c("s1","s1","s1",NA), ed=c(2,3,4,5), total=rep(10,4),
    conditions=c("A","A","B","B")
  )
  expect_error(
    fit_bb_repeated(dat, site="site_id"),
    "Missing values"
  )
})

test_that("first observed condition is the default reference", {
  dat <- data.frame(
    site_id=rep("s1",4), ed=c(2,3,4,5), total=rep(10,4),
    conditions=c("B","A","B","A")
  )
  fit <- fit_bb_repeated(
    dat, site="site_id", family_order="binomial",
    method="coefficient", min_samples=2, quiet=TRUE
  )
  expect_equal(fit$settings$ref_level, NULL)
  expect_equal(levels(fit$data$.condition__), c("B","A"))
})

test_that("min_samples and min_subjects are combined across conditions", {
  # s1: A=3, B=1; s2: A=1, B=3. Both sites have 4 observations
  # and 2+ unique subjects overall, but neither condition reaches 4.
  dat <- data.frame(
    site_id=rep(c("s1","s2"), each=4),
    ed=c(2,3,4,5, 2,3,4,5),
    total=rep(10,8),
    conditions=c("A","A","A","B", "A","B","B","B"),
    indID=c("S1","S1","S2","S2", "S1","S2","S2","S3")
  )
  fit <- fit_bb_repeated(
    dat, site="site_id", subject="indID",
    family_order="binomial", min_samples=4, min_subjects=2,
    quiet=TRUE
  )
  expect_equal(length(fit$fits), 2)
  expect_true(all(vapply(fit$fits, function(x) !is.null(x$fit), logical(1))))
})


test_that("zero-total observations are rejected through min_total", {
  dat <- data.frame(
    site_id=rep("s1",4), ed=c(0,1,2,3), total=c(0,10,10,10),
    conditions=c("A","A","B","B")
  )
  expect_error(
    fit_bb_repeated(dat, site="site_id", min_total=0),
    "greater than or equal to 1"
  )
})

test_that("invalid contrast levels are handled per site", {
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
  expect_true(!is.na(fit$fits[["s1"]]$engine))
  expect_equal(fit$fits[["s2"]]$status, "not_contrastable")
  expect_true(is.null(fit$fits[["s2"]]$fit))
  expect_false(fit$fits[["s2"]]$fit_attempted)
  expect_match(fit$fits[["s2"]]$failure_reason, "not present")
})



test_that("site-specific condition levels are dropped before fitting", {
  dat <- data.frame(
    site_id=rep(c("s1","s2"), each=4),
    ed=c(2,3,4,5, 2,3,4,5), total=rep(10,8),
    conditions=c("A","A","B","B", "B","B","C","C")
  )
  fit <- fit_bb_repeated(
    dat, site="site_id", ref_level="A",
    family_order="binomial", method="coefficient",
    contrast_pairs=data.frame(A="B", B="C"), quiet=TRUE
  )

  expect_equal(fit$fits[["s1"]]$status, "not_contrastable")
  expect_false(fit$fits[["s1"]]$fit_attempted)
  expect_true(is.null(fit$fits[["s1"]]$fit))

  x <- fit$fits[["s2"]]
  expect_equal(x$status, "success")
  expect_true(x$fit_attempted)
  expect_false(is.null(x$fit))
  expect_equal(levels(model.frame(x$fit)[[".condition__"]]), c("B", "C"))
  expect_equal(x$contrasts$contrast, "B_vs_C")
})

test_that("a site observed in only one condition is retained with a structured failure result", {
  dat <- data.frame(
    site_id=rep(c("s1","s2"), each=3),
    ed=c(2,3,4, 1,2,3), total=rep(10,6),
    conditions=c("A","A","A", "A","B","B")
  )
  fit <- fit_bb_repeated(
    dat, site="site_id", family_order="binomial",
    min_samples=3, quiet=TRUE
  )
  expect_equal(length(fit$fits), 2)
  expect_true(is.null(fit$fits[["s1"]]$fit))
  expect_equal(fit$fits[["s1"]]$engine, NA_character_)
  expect_equal(fit$fits[["s1"]]$status, "not_contrastable")
  expect_equal(fit$fits[["s1"]]$failure_reason, "fewer than two condition levels")
  expect_equal(nrow(fit$fits[["s1"]]$contrasts), 0L)
})


test_that("binomial-first selection is not labeled as fallback and diagnostics are retained", {
  dat <- data.frame(
    site_id=rep("s1", 4), ed=c(2,3,4,5), total=rep(10,4),
    conditions=c("A","A","B","B")
  )
  fit <- fit_bb_repeated(
    dat, site="site_id", family_order="binomial",
    method="coefficient", quiet=TRUE
  )
  x <- fit$fits[[1]]
  expect_equal(x$engine, "binomial")
  expect_false(x$fallback)
  expect_true("binomial" %in% names(x$diagnostics_by_family))
  expect_equal(x$diagnostics_by_family$binomial$fit_error, NULL)
})

test_that("prediction uses the fitted model", {
  dat <- data.frame(
    site_id=rep("s1",4), ed=c(2,3,4,5), total=rep(10,4),
    conditions=c("A","A","B","B")
  )
  fit <- fit_bb_repeated(
    dat, site="site_id", family_order="binomial",
    method="coefficient", quiet=TRUE
  )
  if (!is.null(fit$fits[[1]]$fit)) {
    pr <- predict_bb_repeated(fit$fits[[1]]$fit)
    expect_equal(nrow(pr), 2L)
    expect_true(all(is.finite(pr$predicted_probability)))
  }
})


test_that("factor count columns are rejected rather than coerced to factor codes", {
  dat <- data.frame(
    site_id=rep("s1",4), ed=factor(c("2","3","4","5")), total=rep(10,4),
    conditions=c("A","A","B","B")
  )
  expect_error(
    fit_bb_repeated(dat, site="site_id", family_order="binomial", quiet=TRUE),
    "must be numeric count columns"
  )
})


test_that("fractional counts are rejected", {
  dat <- data.frame(
    site_id=rep("s1",4), ed=c(1.5,2.5,3,4), total=rep(10,4),
    conditions=c("A","A","B","B")
  )
  expect_error(
    fit_bb_repeated(dat, site="site_id", family_order="binomial", quiet=TRUE),
    "integer-valued counts"
  )
})

test_that("beta-binomial fitting path is exercised", {
  skip_if_not_installed("glmmTMB")
  set.seed(123)
  n_subject <- 12
  dat <- data.frame(
    site_id=rep("s1", 40),
    conditions=rep(c("A","B"), each=20),
    total=rep(100L, 40)
  )
  p <- ifelse(dat$conditions=="A", 0.10, 0.20)
  dat$ed <- rbinom(nrow(dat), size=dat$total, prob=p)

  fit <- fit_bb_repeated(
    dat, site="site_id",
    family_order="beta-binomial",
    method="both", min_samples=2,
    p_adjust_method="BH", quiet=TRUE
  )
  expect_length(fit$fits, 1)
  expect_true(!is.null(fit$fits[[1]]$fit))
  expect_equal(fit$fits[[1]]$engine, "beta-binomial")
  expect_true("coefficient_p_adj" %in% names(fit$fits[[1]]$contrasts))
  expect_true("delta_p_adj" %in% names(fit$fits[[1]]$contrasts))
})

test_that("multiple-testing correction is applied separately within each contrast", {
  skip_if_not_installed("glmmTMB")
  set.seed(456)
  dat <- expand.grid(
    site_id=paste0("s",1:3),
    conditions=c("A","B","C"),
    indID=paste0("S",1:8)
  )
  dat$total <- 30L
  dat$ed <- rbinom(
    nrow(dat), dat$total,
    ifelse(dat$conditions=="A", 0.10,
           ifelse(dat$conditions=="B", 0.15, 0.30))
  )
  fit <- fit_bb_repeated(
    dat, site="site_id", subject="indID",
    family_order="binomial", method="both",
    min_samples=2, min_subjects=2, p_adjust_method="BH", quiet=TRUE
  )

  all_contrasts <- do.call(rbind, lapply(fit$fits, function(x) {
    if (identical(x$status, "success")) x$contrasts else NULL
  }))
  expect_true(all(c("A_vs_B", "A_vs_C", "B_vs_C") %in% all_contrasts$contrast))

  for (ct in unique(all_contrasts$contrast)) {
    z <- all_contrasts[all_contrasts$contrast == ct, , drop=FALSE]
    expect_equal(
      z$coefficient_p_adj,
      p.adjust(z$coefficient_p, method="BH")
    )
    expect_equal(
      z$delta_p_adj,
      p.adjust(z$delta_p, method="BH")
    )
  }
})


test_that("marginal prediction is available as an optional reporting scale", {
  skip_if_not_installed("glmmTMB")
  set.seed(789)
  dat <- data.frame(
    site_id=rep("s1", 40),
    conditions=rep(c("A","A","B","B"), 10),
    indID=rep(paste0("S", 1:10), each=4),
    total=rep(100L, 40)
  )
  p <- ifelse(dat$conditions=="A", 0.20, 0.60)
  dat$ed <- rbinom(nrow(dat), dat$total, p)
  fit <- fit_bb_repeated(
    dat, site="site_id", subject="indID", family_order="binomial",
    method="coefficient", min_samples=2, min_subjects=2, quiet=TRUE
  )
  if (!is.null(fit$fits[[1]]$fit)) {
    pr <- predict_bb_repeated(fit$fits[[1]]$fit, type="both")
    expect_true(all(c("conditional_probability", "marginal_probability") %in% names(pr)))
    expect_equal(nrow(pr), 2L)
    expect_true(all(is.finite(pr$conditional_probability)))
    expect_true(all(is.finite(pr$marginal_probability)))
    expect_true(all(pr$marginal_probability >= 0 & pr$marginal_probability <= 1))
  }
})
