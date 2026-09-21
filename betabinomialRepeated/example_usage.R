
set.seed(2026)
dat <- expand.grid(
  site_id = paste0("site", 1:10),
  indID = paste0("ID", 1:20),
  conditions = c("HL","ND","RM","PO")
)
dat$total <- sample(10:40, nrow(dat), replace=TRUE)
prob <- c(HL=.05, ND=.08, RM=.12, PO=.18)
dat$ed <- rbinom(nrow(dat), dat$total, prob[dat$conditions])

fit <- fit_bb_repeated(
  dat, site="site_id", edited="ed", total="total",
  condition="conditions", subject="indID",
  ref_level="HL", method="both"
)

# All site/contrast results:
all_results <- dplyr::bind_rows(
  lapply(names(fit$fits), function(s) {
    z <- fit$fits[[s]]
    if (is.null(z) || !nrow(z$contrasts)) return(NULL)
    dplyr::mutate(z$contrasts, site_id=s)
  })
)

comparison <- dplyr::bind_rows(
  lapply(names(fit$fits), function(s) {
    z <- fit$fits[[s]]
    if (is.null(z) || !nrow(z$comparison)) return(NULL)
    dplyr::mutate(z$comparison, site_id=s)
  })
)

# Optional: require at least 2 independent subjects per site when repeated
# measures are present.
# fit <- fit_bb_repeated(dat, site="site_id", edited="ed", total="total",
#                         condition="conditions", subject="indID",
#                         min_samples=2, min_subjects=2)

# Optional: compare conditional (b = 0) and marginal population-averaged predictions
# predict_bb_repeated(fit$fits[["chr1:123:+"]]$fit, type = "both")
