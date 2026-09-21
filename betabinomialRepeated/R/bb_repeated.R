#' Fit beta-binomial regression with repeated measures and binomial fallback
#'
#' @description
#' Fits one model per site using a beta-binomial mixed model when possible,
#' with automatic binomial fallback. Repeated measures are represented by a
#' random intercept for `subject`.
#'
#' `method = "coefficient"` returns contrasts on the model's log-odds scale.
#' `method = "delta"` returns contrasts of fitted probabilities using the
#' delta method. When a subject random intercept is present, these are
#' conditional probabilities evaluated at a random effect of zero, not
#' population-averaged (marginal) probabilities. `method = "both"` returns
#' both scales and a direct comparison.
#'
#' By default, raw Wald p-values are adjusted using the Benjamini-Hochberg
#' (BH) procedure separately within each pairwise contrast and separately for
#' coefficient- and delta-method p-values.
#'
#' @param data A data.frame/tibble with one row per site-sample observation.
#' @param site Column name identifying sites. If NULL, all rows are treated as
#' one site.
#' @param edited Column name containing edited-read counts.
#' @param total Column name containing total-read counts.
#' @param condition Column name containing the experimental condition.
#' @param subject Optional column name identifying repeated-measure subjects.
#' @param ref_level Reference condition. If NULL, the first observed condition
#' is used as the reference level.
#' @param method One of `"coefficient"`, `"delta"`, or `"both"`.
#' @param family_order Families attempted in the specified order. The first
#' family is the primary method; subsequent families are fallback methods if
#' an earlier family is not accepted.
#' @param min_total Minimum positive total reads required for a row to enter a
#' model. Must be greater than or equal to 1.
#' @param min_samples Minimum number of qualifying observations required across
#' all conditions combined at a site. Filtering is not condition-specific, and
#' a site may be retained even if it is observed in only one condition.
#' @param min_subjects Optional minimum number of unique subjects required across
#' all conditions combined at a site. Ignored when `subject = NULL`.
#' @param conf_level Confidence level for intervals.
#' @param contrast_pairs Optional two-column data.frame/matrix defining
#' comparisons. By default all unordered pairs are generated and reported as
#' `A_vs_B` with the first group minus the second group.
#' @param p_adjust_method Multiple-testing adjustment method passed to
#' [stats::p.adjust()]. Default is `"BH"`. Use `"none"` to retain raw p-values.
#' Adjustment is performed separately within each pairwise contrast and
#' separately for coefficient- and delta-method p-values.
#' @param control A `glmmTMBControl()` object controlling model optimization.
#' @param quiet Suppress per-site messages.
#' @return An object of class `bb_repeated`.
#' @examples
#' \dontrun{
#' fit <- fit_bb_repeated(dat, site="site_id", edited="ed", total="total",
#'                        condition="conditions", subject="indID",
#'                        ref_level="HL", method="both")
#' fit$fits[["chr1:123:+"]]$contrasts
#' fit$fits[["chr1:123:+"]]$comparison
#' }
#' @export
fit_bb_repeated <- function(data, site = NULL, edited = "ed", total = "total",
                            condition = "conditions", subject = NULL,
                            ref_level = NULL, method = c("both","coefficient","delta"),
                            family_order = c("beta-binomial","binomial"),
                            min_total = 1L, min_samples = 2L,
                            min_subjects = NULL,
                            conf_level = 0.95, contrast_pairs = NULL,
                            p_adjust_method = "BH",
                            control = glmmTMB::glmmTMBControl(), quiet = FALSE) {
  method <- match.arg(method)
  stopifnot(is.data.frame(data))
  if (!is.numeric(min_total) || length(min_total) != 1L || !is.finite(min_total) || min_total < 1)
    stop("min_total must be a single number greater than or equal to 1.")
  if (!is.numeric(min_samples) || length(min_samples) != 1L || !is.finite(min_samples) || min_samples < 1)
    stop("min_samples must be a single positive number.")
  if (!is.null(min_subjects) &&
      (!is.numeric(min_subjects) || length(min_subjects) != 1L ||
       !is.finite(min_subjects) || min_subjects < 1))
    stop("min_subjects must be NULL or a single positive number.")
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be between 0 and 1.")
  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      !(p_adjust_method %in% stats::p.adjust.methods))
    stop("p_adjust_method must be one of stats::p.adjust.methods.")
  if (!all(family_order %in% c("beta-binomial", "binomial")) ||
      length(family_order) < 1L)
    stop("family_order must contain only 'beta-binomial' and/or 'binomial'.")
  family_order <- unique(family_order)

  req <- c(edited, total, condition)
  if (!is.null(site)) req <- c(req, site)
  if (!is.null(subject)) req <- c(req, subject)
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "))

  d <- tibble::as_tibble(data)

  # Counts must be numeric, finite, non-negative, and integer-valued.
  # In particular, do not coerce factors with as.numeric(), because that
  # returns the underlying factor codes rather than the displayed counts.
  if (!is.numeric(d[[edited]]) || !is.numeric(d[[total]]))
    stop("edited and total must be numeric count columns.")
  d$.edited__ <- d[[edited]]
  d$.total__ <- d[[total]]

  cat_cols <- c(condition, site, subject)
  cat_cols <- cat_cols[!is.na(cat_cols)]
  if (any(vapply(cat_cols, function(x) any(is.na(d[[x]])), logical(1)))) {
    bad <- cat_cols[vapply(cat_cols, function(x) any(is.na(d[[x]])), logical(1))]
    stop("Missing values are not allowed in categorical identifier columns: ",
         paste(bad, collapse=", "))
  }

  if (any(!is.finite(d$.edited__)) || any(!is.finite(d$.total__)))
    stop("Counts must be finite.")
  if (any(d$.edited__ != floor(d$.edited__)) ||
      any(d$.total__ != floor(d$.total__)))
    stop("edited and total must contain integer-valued counts.")
  if (any(d$.edited__ < 0 | d$.total__ < 0 | d$.edited__ > d$.total__))
    stop("Require 0 <= edited <= total for every row.")

  # Preserve the first observed condition as the default reference level.
  condition_levels <- unique(as.character(d[[condition]]))
  d$.condition__ <- factor(d[[condition]], levels=condition_levels)
  if (!is.null(ref_level)) {
    if (!ref_level %in% levels(d$.condition__))
      stop("ref_level is not a condition level.")
    d$.condition__ <- stats::relevel(d$.condition__, ref = ref_level)
  }
  if (!is.null(subject)) d$.subject__ <- factor(d[[subject]])

  if (is.null(site)) {
    d$.site__ <- factor("all_sites")
  } else {
    d$.site__ <- factor(d[[site]])
  }

  # Filtering is performed across all conditions combined within each site.
  # Sites are retained based on total qualifying observations/subjects, even
  # when a site is observed in only one condition.
  keep <- d |>
    dplyr::filter(.total__ >= min_total) |>
    dplyr::group_by(.site__) |>
    dplyr::filter(
      dplyr::n() >= min_samples,
      if (is.null(min_subjects) || is.null(subject)) TRUE
      else dplyr::n_distinct(.subject__) >= min_subjects
    ) |>
    dplyr::ungroup()

  sites <- levels(droplevels(keep$.site__))
  fits <- setNames(vector("list", length(sites)), sites)
  for (s in sites) {
    if (!quiet) message("Fitting site: ", s)
    dd <- dplyr::filter(keep, .site__ == s)
    # A site may contain only a subset of the global condition levels.
    # Drop unused levels before site-wise model fitting so that the model
    # matrix reflects the conditions actually observed at this site.
    dd$.condition__ <- droplevels(dd$.condition__)
    fits[[s]] <- fit_bb_repeated_site(
      dd, condition_col=".condition__", edited_col=".edited__",
      total_col=".total__", subject_col=if (is.null(subject)) NULL else ".subject__",
      method=method, family_order=family_order, conf_level=conf_level,
      contrast_pairs=contrast_pairs, control=control, quiet=quiet
    )
  }

  # Apply multiplicity correction separately within each contrast and
  # separately for coefficient- and delta-method p-values. Thus, for example,
  # all site-level p-values for A_vs_B form one BH family, while A_vs_C forms
  # a separate family. This matches the intended analysis, where each
  # pairwise contrast is treated as its own multiple-testing family.
  if (length(fits)) {
    contrast_names <- unique(unlist(lapply(fits, function(x) {
      if (is.null(x) || is.null(x$contrasts) || !nrow(x$contrasts))
        return(character())
      as.character(x$contrasts$contrast)
    })))

    for (ct in contrast_names) {
      if (method %in% c("coefficient", "both")) {
        vals <- lapply(fits, function(x) {
          if (is.null(x) || is.null(x$contrasts) ||
              !"coefficient_p" %in% names(x$contrasts))
            return(numeric())
          x$contrasts$coefficient_p[x$contrasts$contrast == ct]
        })
        p <- unlist(vals, use.names=FALSE)
        if (length(p)) {
          adj <- stats::p.adjust(p, method=p_adjust_method)
          k <- 1L
          for (i in seq_along(fits)) {
            x <- fits[[i]]
            if (is.null(x) || is.null(x$contrasts) ||
                !"coefficient_p" %in% names(x$contrasts))
              next
            idx <- which(x$contrasts$contrast == ct)
            if (!length(idx))
              next
            x$contrasts$coefficient_p_adj[idx] <- adj[k:(k+length(idx)-1L)]
            fits[[i]] <- x
            k <- k + length(idx)
          }
        }
      }

      if (method %in% c("delta", "both")) {
        vals <- lapply(fits, function(x) {
          if (is.null(x) || is.null(x$contrasts) ||
              !"delta_p" %in% names(x$contrasts))
            return(numeric())
          x$contrasts$delta_p[x$contrasts$contrast == ct]
        })
        p <- unlist(vals, use.names=FALSE)
        if (length(p)) {
          adj <- stats::p.adjust(p, method=p_adjust_method)
          k <- 1L
          for (i in seq_along(fits)) {
            x <- fits[[i]]
            if (is.null(x) || is.null(x$contrasts) ||
                !"delta_p" %in% names(x$contrasts))
              next
            idx <- which(x$contrasts$contrast == ct)
            if (!length(idx))
              next
            x$contrasts$delta_p_adj[idx] <- adj[k:(k+length(idx)-1L)]
            fits[[i]] <- x
            k <- k + length(idx)
          }
        }
      }
    }

    for (i in seq_along(fits)) {
      x <- fits[[i]]
      if (!is.null(x) && !is.null(x$contrasts) && nrow(x$contrasts) &&
          method == "both")
        x$comparison <- compare_bb_methods(x$contrasts)
      fits[[i]] <- x
    }
  }

  site_status <- vapply(fits, function(x) {
    if (is.null(x) || is.null(x$status)) "model_failed" else x$status
  }, character(1))
  n_sites_retained <- length(fits)
  n_sites_attempted <- sum(vapply(fits, function(x) isTRUE(x$fit_attempted), logical(1)))
  n_sites_successful <- sum(site_status == "success")
  n_sites_failed <- sum(site_status == "model_failed")
  n_sites_not_contrastable <- sum(site_status == "not_contrastable")

  out <- list(
    call = match.call(), fits = fits, data = keep,
    n_sites_retained = n_sites_retained,
    n_sites_attempted = n_sites_attempted,
    n_sites_successful = n_sites_successful,
    n_sites_failed = n_sites_failed,
    n_sites_not_contrastable = n_sites_not_contrastable,
    settings = list(site=site, edited=edited, total=total,
                    condition=condition, subject=subject, ref_level=ref_level,
                    method=method, family_order=family_order,
                    min_total=min_total, min_samples=min_samples,
                    min_subjects=min_subjects, conf_level=conf_level,
                    p_adjust_method=p_adjust_method, control=control)
  )
  class(out) <- "bb_repeated"
  out
}

#' @rdname fit_bb_repeated
#' @export
fit_bb_repeated_site <- function(data, condition_col=".condition__",
                                 edited_col=".edited__", total_col=".total__",
                                 subject_col=NULL,
                                 method=c("both","coefficient","delta"),
                                 family_order=c("beta-binomial","binomial"),
                                 conf_level=0.95, contrast_pairs=NULL,
                                 control=glmmTMB::glmmTMBControl(), quiet=TRUE) {
  method <- match.arg(method)
  n_conditions <- dplyr::n_distinct(data[[condition_col]])
  n_subjects <- if (is.null(subject_col)) NA_integer_ else dplyr::n_distinct(data[[subject_col]])
  n_observations <- nrow(data)

  # Site-level data can inherit factor levels from the full analysis even
  # when some conditions are not observed at this site. Drop those unused
  # levels before constructing the model, so the model matrix contains only
  # estimable site-level condition effects.
  data[[condition_col]] <- droplevels(factor(data[[condition_col]]))
  levs <- levels(data[[condition_col]])
  n_conditions <- length(levs)

  if (n_conditions < 2L) {
    return(list(
      status="not_contrastable",
      fit_attempted=FALSE,
      engine=NA_character_, fit=NULL, errors=character(),
      diagnostics=NULL, diagnostics_by_family=list(),
      fallback=FALSE, fallback_reason=NA_character_,
      failure_reason="fewer than two condition levels",
      n_observations=n_observations, n_subjects=n_subjects,
      n_conditions=n_conditions,
      contrasts=tibble::tibble(), comparison=tibble::tibble()
    ))
  }

  # Validate requested contrasts before attempting a model fit. A site that
  # lacks a requested condition is structurally not contrastable, and should
  # not incur an optimizer attempt merely to discover that fact.
  requested_contrasts <- contrast_pairs
  if (is.null(requested_contrasts)) {
    cc <- utils::combn(levs, 2, simplify=FALSE)
    requested_contrasts <- do.call(rbind, lapply(cc, function(x) c(x[1], x[2])))
  }
  requested_contrasts <- as.data.frame(requested_contrasts, stringsAsFactors=FALSE)
  if (ncol(requested_contrasts) != 2)
    stop("contrast_pairs must have exactly two columns.")
  names(requested_contrasts) <- c("A","B")

  invalid <- is.na(requested_contrasts$A) | is.na(requested_contrasts$B) |
             !requested_contrasts$A %in% levs | !requested_contrasts$B %in% levs
  if (any(invalid)) {
    missing_levels <- unique(c(
      as.character(requested_contrasts$A[is.na(requested_contrasts$A) |
                                         !requested_contrasts$A %in% levs]),
      as.character(requested_contrasts$B[is.na(requested_contrasts$B) |
                                         !requested_contrasts$B %in% levs])
    ))
    missing_levels <- missing_levels[!is.na(missing_levels)]
    reason <- if (length(missing_levels)) {
      paste0("requested contrast levels not present in site: ",
             paste(missing_levels, collapse=", "))
    } else {
      "requested contrast contains missing condition levels"
    }
    return(list(
      status="not_contrastable",
      fit_attempted=FALSE,
      engine=NA_character_, fit=NULL, errors=character(),
      diagnostics=NULL, diagnostics_by_family=list(),
      fallback=FALSE, fallback_reason=NA_character_,
      failure_reason=reason,
      n_observations=n_observations, n_subjects=n_subjects,
      n_conditions=n_conditions,
      contrasts=tibble::tibble(), comparison=tibble::tibble()
    ))
  }

  form <- if (is.null(subject_col))
    stats::as.formula(sprintf("cbind(%s, %s-%s) ~ %s", edited_col, total_col,
                              edited_col, condition_col))
  else
    stats::as.formula(sprintf("cbind(%s, %s-%s) ~ %s + (1|%s)",
                              edited_col, total_col, edited_col, condition_col, subject_col))

  fit <- NULL
  engine <- NULL
  errors <- character()
  diagnostics <- NULL
  diagnostics_by_family <- list()
  fallback <- FALSE
  fallback_reason <- NA_character_

  for (fam in unique(family_order)) {
    f <- if (fam == "beta-binomial") glmmTMB::betabinomial() else stats::binomial()
    z <- try(glmmTMB::glmmTMB(formula=form, family=f, data=data, control=control), silent=TRUE)

    if (inherits(z, "try-error")) {
      err_msg <- as.character(z)
      errors <- c(errors, paste0(fam, ": ", err_msg))
      diagnostics_by_family[[fam]] <- list(
        convergence=NA_integer_, pdHess=NA,
        convergence_ok=FALSE, hessian_ok=FALSE,
        fit_error=err_msg
      )
      next
    }

    conv <- z$fit$convergence
    pd_hess <- if (!is.null(z$sdr$pdHess)) isTRUE(z$sdr$pdHess) else NA
    fam_diag <- list(
      convergence=conv,
      pdHess=pd_hess,
      convergence_ok=identical(as.integer(conv), 0L),
      hessian_ok=isTRUE(pd_hess),
      fit_error=NULL
    )
    diagnostics_by_family[[fam]] <- fam_diag

    if (!fam_diag$convergence_ok || !fam_diag$hessian_ok) {
      reason_parts <- c(
        if (!fam_diag$convergence_ok) "nonzero optimizer convergence" else NULL,
        if (!fam_diag$hessian_ok) "non-positive-definite Hessian" else NULL
      )
      reason <- paste(reason_parts, collapse="; ")
      errors <- c(errors, paste0(fam, ": ", reason))
      diagnostics_by_family[[fam]]$fit_error <- reason
      next
    }

    fit <- z
    diagnostics <- fam_diag
    engine <- fam
    break
  }

  if (is.null(fit)) {
    return(list(
      status="model_failed",
      fit_attempted=TRUE,
      engine=NA_character_, fit=NULL, errors=errors,
      diagnostics=diagnostics, diagnostics_by_family=diagnostics_by_family,
      fallback=FALSE, fallback_reason=NA_character_,
      failure_reason="all requested model families failed",
      n_observations=n_observations, n_subjects=n_subjects,
      n_conditions=n_conditions,
      contrasts=tibble::tibble(), comparison=tibble::tibble()
    ))
  }

  # The first requested family is primary; later families are fallbacks.
  if (!identical(engine, family_order[[1L]])) {
    fallback <- TRUE
    first_family <- family_order[[1L]]
    first_diag <- diagnostics_by_family[[first_family]]
    fallback_reason <- if (!is.null(first_diag$fit_error)) {
      first_diag$fit_error
    } else {
      reason_parts <- c(
        if (!isTRUE(first_diag$convergence_ok)) "nonzero optimizer convergence" else NULL,
        if (!isTRUE(first_diag$hessian_ok)) "non-positive-definite Hessian" else NULL
      )
      if (length(reason_parts)) paste(reason_parts, collapse="; ")
      else paste0(first_family, " was not accepted.")
    }
  }

  contrast_pairs <- requested_contrasts

  res <- extract_bb_contrasts(
    fit, contrast_pairs, method=method, conf_level=conf_level, engine=engine
  )
  comp <- if (method == "both") compare_bb_methods(res) else tibble::tibble()

  list(
    status="success",
    fit_attempted=TRUE,
    engine=engine, fit=fit, diagnostics=diagnostics,
    diagnostics_by_family=diagnostics_by_family,
    fallback=fallback, fallback_reason=fallback_reason,
    failure_reason=NA_character_,
    n_observations=n_observations, n_subjects=n_subjects,
    n_conditions=n_conditions,
    contrasts=res, comparison=comp
  )
}

#' Extract coefficient- and/or delta-method contrasts
#'
#' @param fit A glmmTMB fit.
#' @param contrast_pairs Two-column data.frame with A and B condition names.
#' @param conf_level Confidence level.
#' @param engine Model engine label.
#' @return A tibble containing requested contrast scales. Delta-method
#' probabilities are conditional probabilities evaluated at a random effect of zero
#' when a subject random intercept is present; they are not
#' population-averaged probabilities. Confidence intervals for delta-method
#' probability differences are Wald intervals and are not constrained to the
#' mathematical range [-1, 1].
#' @export
extract_bb_contrasts <- function(fit, contrast_pairs,
                                 method=c("both","coefficient","delta"),
                                 conf_level=.95, engine=NULL) {
  method <- match.arg(method)
  mf <- stats::model.frame(fit)
  cond_terms <- attr(stats::terms(fit), "term.labels")
  cond_terms <- cond_terms[!grepl("^\\(1\\|", cond_terms)]
  cond_terms <- cond_terms[!grepl("\\|", cond_terms)]
  if (length(cond_terms) != 1L)
    stop("Could not uniquely identify the condition term in the fitted model.")
  condition_var <- cond_terms[[1L]]

  beta <- glmmTMB::fixef(fit)$cond
  if (is.null(beta) || !length(beta) || is.null(names(beta))) {
    stop("Could not extract conditional fixed effects from the fitted model.")
  }

  ## Extract the conditional fixed-effect covariance matrix and explicitly
  ## align it to the fixed-effect coefficient names. This is important for
  ## mixed models because vcov() may return multiple parameter blocks.
  V <- stats::vcov(fit)
  if (is.list(V)) {
    if (is.null(V$cond)) {
      stop("Could not extract the conditional fixed-effect covariance matrix.")
    }
    V <- V$cond
  }
  if (!is.matrix(V) || is.null(rownames(V)) || is.null(colnames(V))) {
    stop("Conditional fixed-effect covariance matrix is not a named matrix.")
  }
  if (!all(names(beta) %in% rownames(V)) ||
      !all(names(beta) %in% colnames(V))) {
    stop("Conditional fixed-effect covariance matrix does not contain all fixed-effect coefficients.")
  }
  V <- V[names(beta), names(beta), drop=FALSE]

  cond_factor <- mf[[condition_var]]
  if (is.null(cond_factor) || !is.factor(cond_factor))
    stop("The fitted condition variable must be a factor.")
  levs <- levels(cond_factor)

  if (any(!contrast_pairs$A %in% levs) || any(!contrast_pairs$B %in% levs))
    stop("contrast_pairs contains levels not present in the fitted model.")

  make_newdata <- function(level) {
    nd <- mf[1, , drop=FALSE]
    nd[[condition_var]] <- factor(level, levels=levs)
    nd
  }
  ## The regression remains the full repeated-measures model, including
  ## (1 | subject). For delta-method calculations, however, only the
  ## conditional fixed-effects design matrix is required because beta is
  ## extracted from fixef(fit)$cond. Do not reconstruct X from the full
  ## mixed-model terms object.
  fixed_formula <- stats::reformulate(condition_var)
  fixed_contrasts <- attr(mf, "contrasts")
  x_for <- function(level) {
    nd <- make_newdata(level)
    X <- stats::model.matrix(
      fixed_formula,
      nd,
      contrasts.arg=fixed_contrasts
    )
    if (!all(names(beta) %in% colnames(X))) {
      stop("Fixed-effect design matrix does not contain all fitted coefficients.")
    }
    X <- X[, names(beta), drop=FALSE]
    as.numeric(X[1, ])
  }
  eta <- function(level) {
    x <- x_for(level)
    e <- sum(x * beta)
    p <- stats::plogis(e)
    list(x=x, eta=e, p=p)
  }
  zcrit <- stats::qnorm(1-(1-conf_level)/2)

  rows <- lapply(seq_len(nrow(contrast_pairs)), function(i) {
    A <- as.character(contrast_pairs$A[i]); B <- as.character(contrast_pairs$B[i])
    eA <- eta(A); eB <- eta(B)
    out <- tibble::tibble(A=A, B=B, contrast=paste0(A,"_vs_",B), engine=engine)

    if (method %in% c("delta","both")) {
      gA <- eA$p*(1-eA$p)*eA$x
      gB <- eB$p*(1-eB$p)*eB$x
      gd <- gA-gB
      se_delta <- sqrt(max(0, as.numeric(t(gd) %*% V %*% gd)))
      pd <- eA$p-eB$p
      z_delta <- if (se_delta > 0) pd/se_delta else NA_real_
      p_delta <- if (is.finite(z_delta)) 2*stats::pnorm(-abs(z_delta)) else NA_real_
      out <- dplyr::bind_cols(out, tibble::tibble(
        A_pred=eA$p, B_pred=eB$p, delta_diff=pd, delta_SE=se_delta,
        delta_z=z_delta, delta_p=p_delta,
        delta_ci_low=pd-zcrit*se_delta, delta_ci_high=pd+zcrit*se_delta
      ))
    }

    if (method %in% c("coefficient","both")) {
      gc <- eA$x-eB$x
      coef_diff <- sum(gc*beta)
      se_coef <- sqrt(max(0, as.numeric(t(gc) %*% V %*% gc)))
      z_coef <- if (se_coef > 0) coef_diff/se_coef else NA_real_
      p_coef <- if (is.finite(z_coef)) 2*stats::pnorm(-abs(z_coef)) else NA_real_
      out <- dplyr::bind_cols(out, tibble::tibble(
        coefficient_diff=coef_diff, coefficient_SE=se_coef,
        coefficient_z=z_coef, coefficient_p=p_coef,
        coefficient_ci_low=coef_diff-zcrit*se_coef,
        coefficient_ci_high=coef_diff+zcrit*se_coef
      ))
    }
    out
  })
  dplyr::bind_rows(rows)
}

#' Compare coefficient and delta-method results
#'
#' @param contrasts Output from [extract_bb_contrasts()].
#' @return A tibble with method differences and explicit raw/adjusted
#' significance indicators. `significance_agreement` uses adjusted p-values
#' when they are available; otherwise it uses raw p-values.
#' @export
compare_bb_methods <- function(contrasts) {
  if (!nrow(contrasts)) return(tibble::tibble())

  out <- contrasts |>
    dplyr::mutate(
      abs_difference=abs(coefficient_diff-delta_diff),
      same_direction=sign(coefficient_diff)==sign(delta_diff)
    )

  if ("coefficient_p" %in% names(out))
    out$coefficient_significant_raw <- out$coefficient_p < 0.05
  if ("delta_p" %in% names(out))
    out$delta_significant_raw <- out$delta_p < 0.05
  if ("coefficient_significant_raw" %in% names(out) &&
      "delta_significant_raw" %in% names(out)) {
    out$significance_agreement_raw <-
      out$coefficient_significant_raw == out$delta_significant_raw
  }

  has_coef_adj <- "coefficient_p_adj" %in% names(out)
  has_delta_adj <- "delta_p_adj" %in% names(out)
  if (has_coef_adj)
    out$coefficient_significant_adj <- out$coefficient_p_adj < 0.05
  if (has_delta_adj)
    out$delta_significant_adj <- out$delta_p_adj < 0.05
  if (has_coef_adj && has_delta_adj) {
    out$significance_agreement_adj <-
      out$coefficient_significant_adj == out$delta_significant_adj
    out$significance_agreement <- out$significance_agreement_adj
  } else if ("significance_agreement_raw" %in% names(out)) {
    out$significance_agreement <- out$significance_agreement_raw
  }

  out
}

#' Predict condition-specific editing probabilities
#'
#' @param fit A glmmTMB fit.
#' @param newdata Optional data frame for prediction. If omitted, predictions
#' are generated for every condition level in the fitted model.
#' @param type Prediction scale: `conditional` (fixed-effect probability at
#' random effect zero), `marginal` (integrated over the fitted random-intercept
#' distribution), or `both`.
#' @return Tibble of condition and predicted probabilities.
#' @details
#' The fitted regression model is unchanged. When a subject random intercept
#' is present, it remains part of the glmmTMB model used to estimate the fixed
#' effects and their covariance. Conditional predictions use the fixed-effect
#' linear predictor with the random effect set to zero. Marginal predictions
#' numerically integrate the conditional response probability over the
#' estimated normal random-intercept distribution. Marginal predictions are
#' descriptive population-averaged probabilities; this function does not
#' calculate marginal-effect standard errors or confidence intervals.
#' @export
predict_bb_repeated <- function(fit, newdata=NULL,
                                type=c("conditional", "marginal", "both")) {
  if (!inherits(fit, "glmmTMB")) stop("fit must be a glmmTMB model.")
  type <- match.arg(type)
  mf <- stats::model.frame(fit)
  terms_obj <- stats::terms(fit)
  term_labels <- attr(terms_obj, "term.labels")
  condition_terms <- term_labels[!grepl("\\|", term_labels)]
  if (length(condition_terms) != 1L)
    stop("Could not uniquely identify the condition term in the fitted model.")
  condition_var <- condition_terms[[1L]]
  condition_factor <- mf[[condition_var]]
  if (!is.factor(condition_factor))
    stop("The fitted condition variable must be a factor.")

  if (is.null(newdata)) {
    nd <- mf[rep(1, nlevels(condition_factor)), , drop=FALSE]
    nd[[condition_var]] <- factor(levels(condition_factor),
                                  levels=levels(condition_factor))
  } else {
    nd <- as.data.frame(newdata)
    if (!condition_var %in% names(nd))
      stop("newdata must contain the fitted condition variable.")
    nd[[condition_var]] <- factor(nd[[condition_var]],
                                  levels=levels(condition_factor))
    if (anyNA(nd[[condition_var]]))
      stop("newdata contains condition levels not present in the fitted model.")
  }

  conditional <- as.numeric(stats::predict(
    fit, newdata=nd, type="response", re.form=NA
  ))

  if (type == "conditional") {
    return(tibble::as_tibble(nd) |>
      dplyr::mutate(predicted_probability=conditional))
  }

  vc <- glmmTMB::VarCorr(fit)$cond
  if (is.null(vc) || !length(vc)) {
    marginal <- conditional
  } else {
    sd_candidates <- unlist(lapply(vc, function(x) attr(x, "stddev")),
                            use.names=FALSE)
    sd_candidates <- sd_candidates[is.finite(sd_candidates)]
    if (!length(sd_candidates))
      stop("Could not extract the random-intercept standard deviation.")
    sd_re <- as.numeric(sd_candidates[1L])
    if (sd_re <= 0) {
      marginal <- conditional
    } else {
      p_safe <- pmin(pmax(conditional, .Machine$double.eps),
                     1 - .Machine$double.eps)
      eta <- stats::qlogis(p_safe)
      marginal <- vapply(eta, function(e) {
        stats::integrate(
          function(b) stats::plogis(e + b) * stats::dnorm(b, mean=0, sd=sd_re),
          lower=-Inf, upper=Inf, rel.tol=1e-10
        )$value
      }, numeric(1))
    }
  }

  out <- tibble::as_tibble(nd)
  if (type == "marginal") {
    out$predicted_probability <- marginal
    return(out)
  }
  out$conditional_probability <- conditional
  out$marginal_probability <- marginal
  out
}

#' @export
summary.bb_repeated <- function(object, ...) {
  fits <- object$fits
  site_status <- vapply(fits, function(x) {
    if (is.null(x) || is.null(x$status)) "model_failed" else x$status
  }, character(1))

  site_table <- tibble::tibble(
    site_id=names(fits),
    status=site_status,
    engine=vapply(fits, function(x)
      if (is.null(x)) NA_character_ else x$engine %||% NA_character_,
      character(1)),
    fallback=vapply(fits, function(x)
      if (is.null(x)) NA else isTRUE(x$fallback), logical(1)),
    fallback_reason=vapply(fits, function(x)
      if (is.null(x) || is.null(x$fallback_reason) ||
          length(x$fallback_reason)==0L || is.na(x$fallback_reason))
        NA_character_ else x$fallback_reason, character(1)),
    failure_reason=vapply(fits, function(x)
      if (is.null(x) || is.null(x$failure_reason) ||
          length(x$failure_reason)==0L || is.na(x$failure_reason))
        NA_character_ else x$failure_reason, character(1)),
    convergence=vapply(fits, function(x)
      if (is.null(x) || is.null(x$diagnostics) ||
          is.null(x$diagnostics$convergence)) NA_integer_
      else as.integer(x$diagnostics$convergence), integer(1)),
    pdHess=vapply(fits, function(x)
      if (is.null(x) || is.null(x$diagnostics) ||
          is.null(x$diagnostics$pdHess)) NA
      else isTRUE(x$diagnostics$pdHess), logical(1)),
    n_observations=vapply(fits, function(x)
      if (is.null(x) || is.null(x$n_observations)) NA_integer_
      else as.integer(x$n_observations), integer(1)),
    n_subjects=vapply(fits, function(x)
      if (is.null(x) || is.null(x$n_subjects)) NA_integer_
      else as.integer(x$n_subjects), integer(1)),
    n_conditions=vapply(fits, function(x)
      if (is.null(x) || is.null(x$n_conditions)) NA_integer_
      else as.integer(x$n_conditions), integer(1)),
    n_contrasts=vapply(fits, function(x)
      if (is.null(x) || is.null(x$contrasts)) 0L else nrow(x$contrasts),
      integer(1)),
    failed=vapply(fits, function(x)
      if (is.null(x)) TRUE else identical(x$status, "model_failed"), logical(1))
  )

  attr(site_table, "site_counts") <- c(
    retained=object$n_sites_retained %||% length(fits),
    attempted=object$n_sites_attempted %||% sum(vapply(fits, function(x) isTRUE(x$fit_attempted), logical(1))),
    successful=object$n_sites_successful %||% sum(site_status == "success"),
    failed=object$n_sites_failed %||% sum(site_status == "model_failed"),
    not_contrastable=object$n_sites_not_contrastable %||% sum(site_status == "not_contrastable")
  )
  site_table
}

`%||%` <- function(a,b) if (is.null(a) || length(a)==0 || is.na(a)) b else a
