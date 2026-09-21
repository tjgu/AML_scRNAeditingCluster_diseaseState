#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(glmmTMB)
  library(broom.mixed)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(ggrepel)
})

options(stringsAsFactors = FALSE)

# =========================
# CONFIG
# =========================
matrix_dir   <- "./data"
matrix_file  <- "editing_matrix.txt"
bam_order    <- file.path(matrix_dir, "bamfile_lists.txt")
meta_file    <- file.path(matrix_dir, "meta_samples.txt")

min_cov       <- 5
min_samples   <- 2
out_prefix    <- "./results/de_bulk_sites/de_bulk_sites"
ref_levels    <- c("HL","ND","RM","PO")  # HL=Healthy, ND=Diag, RM=PostT, PO=PoorOutcome

FDR_THRESH    <- 0.10
OBS_DIFF_MIN  <- 0.05
VOLC_Y_CAP    <- 10
EPS_P         <- 1e-10

# =========================
# INPUTS
# =========================
# Ensure output dir exists (right after CONFIG block)
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(meta_file), file.exists(bam_order), file.exists(file.path(matrix_dir, matrix_file)))

meta <- fread(meta_file, sep = "\t", header = TRUE) %>%
  mutate(conditions = factor(conditions, levels = ref_levels))
stopifnot("seqID" %in% names(meta), "conditions" %in% names(meta))

bam <- fread(bam_order, header = FALSE)$V1 %>% trimws()
stopifnot(length(bam) > 0)

sel_idx <- which(bam %in% meta$seqID)
stopifnot(length(sel_idx) > 0)
sel_seq <- bam[sel_idx]

meta_sel <- meta %>%
  filter(seqID %in% sel_seq) %>%
  distinct(seqID, .keep_all = TRUE) %>%
  select(seqID, conditions, any_of(c("indID")))

# =========================
# READ & PARSE MATRIX (wide -> long)
# Each sample block: [editing_type, total_reads, edited_reads, editing_ratio]
# We ignore type & ratio; compute uned = total - ed; keep total
# =========================
fpath <- file.path(matrix_dir, matrix_file)
DT <- fread(fpath, sep = "\t", header = FALSE, data.table = TRUE, na.strings = c("NA","","NaN"))

expected_cols <- 2 + 4*length(bam)
if (ncol(DT) != expected_cols) stop(sprintf("Matrix has %d cols, expected %d", ncol(DT), expected_cols))
setnames(DT, 1:2, c("site_id","IGNORE_col2"))

# After reading the matrix
message("Matrix dimensions: ", nrow(DT), " sites x ", ncol(DT), " columns")
message("Expected samples: ", length(bam))
message("Selected samples: ", length(sel_seq))
message("Samples in metadata: ", nrow(meta_sel))

make_long_for <- function(j, seqID) {
  start <- 3 + (j-1)*4
  cols  <- start:(start+3)
  tmp <- DT[, ..cols]
  setnames(tmp, c("etype_IGN","total","ed","ratio"))
  tmp[, `:=`(
    total = suppressWarnings(as.numeric(total)),
    ed    = suppressWarnings(as.numeric(ed)),
    ratio = suppressWarnings(as.numeric(ratio))
  )]
  for (nm in c("total","ed","ratio")) tmp[is.na(get(nm)), (nm) := 0]
  tmp[total < 0 | ed < 0, `:=`(total = pmax(total, 0), ed = pmax(ed, 0))]
  tmp[, uned := pmax(0, total - ed)]
  tmp[, `:=`(
    total = as.integer(round(total)),
    ed    = as.integer(round(ed)),
    uned  = as.integer(round(uned))
  )]
  tmp[, `:=`(seqID = seqID, site_id = DT$site_id)]
  tmp[, .(site_id, seqID, uned, ed, total)]
}

message("Parsing bulk matrix into long format for ", length(sel_idx), " selected samples ...")
long <- data.table::rbindlist(purrr::map2(sel_idx, sel_seq, make_long_for))

# coverage & presence
long[, present := total > min_cov]

# attach condition
long <- data.table::as.data.table(dplyr::left_join(long, meta_sel, by = "seqID"))
long[, conditions := droplevels(factor(conditions, levels = ref_levels))]

# regression sanity
stopifnot(all(long$uned + long$ed == long$total), all(long$total >= 0), all(is.finite(long$total)))

# Keep sites present in >= min_samples with cov > min_cov
keep_sites <- long %>%
  filter(present) %>%
  distinct(site_id, seqID) %>%
  count(site_id, name = "n_present_samples") %>%
  filter(n_present_samples >= min_samples)

dat <- long %>%
  inner_join(keep_sites %>% select(site_id), by = "site_id") %>%
  filter(total > 0)

# split site_id to coords
tmp <- tstrsplit(dat$site_id, ":", fixed = TRUE)
dat$chr    <- tmp[[1]]
dat$pos    <- suppressWarnings(as.integer(tmp[[2]]))
dat$strand <- tmp[[3]]

message(sprintf("Kept %d sites after filter (>%d cov in ≥%d samples).",
                length(unique(dat$site_id)), min_cov, min_samples))


# ==== SAVE FILTERED SITES FROM DE STEP (exact universe used downstream) ====
# keep_sites currently has: site_id, n_present_samples
# Add metadata of thresholds used, then write.
keep_sites$min_cov      <- min_cov
keep_sites$min_samples  <- min_samples

fwrite(keep_sites,
       sprintf("%s_passingSites.tsv", out_prefix), sep = "\t")

# (Optional, large) also save the filtered long matrix rows for these sites
# so the annotation step can attach gene/region per site and join back into a matrix if desired.
long_filtered <- long[site_id %in% keep_sites$site_id]
fwrite(long_filtered,
       sprintf("%s_filteredLong.tsv.gz", out_prefix), sep = "\t")
message("Saved passing sites (and optional filtered long).")


# =========================
# MODEL: beta-binomial per site; fallback to binomial
# =========================
has_indID <- "indID" %in% names(dat) && {
  ids <- dat$indID
  ids <- ids[!is.na(ids)]
  length(ids) > 0 && anyDuplicated(ids) > 0
}
form_noRE <- as.formula(cbind(ed, uned) ~ conditions)
form_RE   <- as.formula(cbind(ed, uned) ~ conditions + (1|indID))

eta_builder <- function(fit) {
  beta <- glmmTMB::fixef(fit)$cond
  Vfull <- stats::vcov(fit); if (is.list(Vfull)) Vfull <- Vfull$cond
  base_intercept <- "(Intercept)" %in% names(beta)

  x_for <- function(level, ref = ref_levels[1]) {
    x <- setNames(numeric(length(beta)), names(beta))
    if (base_intercept) x["(Intercept)"] <- 1
    if (level != ref) {
      nm <- paste0("conditions", level)
      if (nm %in% names(beta)) x[nm] <- 1
    }
    x
  }
  eta_p <- function(level) {
    x <- x_for(level)
    eta <- as.numeric(crossprod(x, beta))
    p <- plogis(eta)
    p <- min(max(p, EPS_P), 1 - EPS_P)
    list(level = level, x = x, eta = eta, p = p)
  }
  list(beta = beta, V = Vfull, x_for = x_for, eta_p = eta_p)
}

fit_one_site <- function(df) {
  df <- droplevels(df)
  if (dplyr::n_distinct(df$conditions) < 2) {
    message("Skipping site - insufficient conditions")
    return(NULL)
  }

  form <- if (has_indID) form_RE else form_noRE
  engine <- "beta-binomial"
  fit <- try(glmmTMB::glmmTMB(form, family = glmmTMB::betabinomial(), data = df), silent = TRUE)
  if (inherits(fit, "try-error")) {
    engine <- "binomial"
    fit <- try(glmmTMB::glmmTMB(form, family = stats::binomial(), data = df), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
  }
  # In fit_one_site function:
  if (inherits(fit, "glmmTMB") && fit$fit$convergence != 0) {
    warning("Model convergence issue for site: ", unique(df$site_id)[1])
  }

  # predicted p-hat per condition
  eb <- eta_builder(fit)
  levs <- levels(df$conditions)
  preds <- map_dfr(levs, ~{
    ep <- eb$eta_p(.x); tibble(condition = .x, p_hat = ep$p)
  }) %>% mutate(engine = engine)

  # pairwise deltas via delta method - unordered (combn gives (a,b) with a lower-index, b higher-index)
  pair_labels <- combn(levs, 2, simplify = FALSE)
  one_diff <- function(a, b) {
    V <- eb$V; if (is.null(V)) return(NULL)
    epA <- eb$eta_p(a); epB <- eb$eta_p(b)
    gA  <- as.numeric(epA$p * (1 - epA$p)) * epA$x
    gB  <- as.numeric(epB$p * (1 - epB$p)) * epB$x
    g   <- gA - gB
    se  <- sqrt(as.numeric(t(g) %*% V %*% g))
    diff <- epA$p - epB$p
    z <- ifelse(se > 0, diff / se, NA_real_)
    p <- ifelse(is.na(z), NA_real_, 2*pnorm(-abs(z)))
    tibble(a = a, b = b, pred_diff_ab = diff, SE_diff_ab = se, z_diff_ab = z, p_diff_ab = p, engine = engine)
  }
  d_unord <- purrr::map_dfr(pair_labels, ~ one_diff(.x[1], .x[2]))

  # Keep ONLY one direction as A_vs_B where A is the higher-index group (b), B is lower (a).
  # So contrast label is A_vs_B = b_vs_a, and diff = p(b)-p(a).
  if (nrow(d_unord)) {
    diffs <- d_unord %>%
      transmute(
        contrast = paste0(b, "_vs_", a),
        pred_diff = -pred_diff_ab,   # p(b)-p(a)
        SE_diff   =  SE_diff_ab,
        z_diff    = -z_diff_ab,
        p_diff    =  p_diff_ab,
        engine    = engine
      )
  } else {
    diffs <- tibble()
  }

  # observed pooled fractions per condition
  obs <- df %>%
    group_by(conditions) %>%
    summarise(
      obs_sum_ed  = sum(ed),
      obs_sum_tot = sum(total),
      obs_frac    = ifelse(obs_sum_tot > 0, obs_sum_ed / obs_sum_tot, NA_real_),
      .groups = "drop"
    )

  list(preds = preds, diffs = diffs, obs = obs)
}

# fit per site
message("Fitting models…")
by_site <- split(dat, dat$site_id)
fits <- lapply(by_site, fit_one_site)

# After model fitting, add:
message("Model fitting summary:")
message("  Beta-binomial models: ", sum(sapply(fits, function(x) if(is.null(x)) FALSE else x$diffs$engine[1] == "beta-binomial"), na.rm = TRUE))
message("  Binomial fallback: ", sum(sapply(fits, function(x) if(is.null(x)) FALSE else x$diffs$engine[1] == "binomial"), na.rm = TRUE))
message("  Failed models: ", sum(sapply(fits, is.null)))

# collect
diffs <- bind_rows(lapply(names(fits), function(s) {
  x <- fits[[s]]; if (is.null(x) || is.null(x$diffs) || nrow(x$diffs)==0) return(NULL)
  cbind(site_id = s, x$diffs)
}))
preds <- bind_rows(lapply(names(fits), function(s) {
  x <- fits[[s]]; if (is.null(x) || is.null(x$preds) || nrow(x$preds)==0) return(NULL)
  cbind(site_id = s, x$preds)
}))
obs_summ <- bind_rows(lapply(names(fits), function(s) {
  x <- fits[[s]]; if (is.null(x) || is.null(x$obs) || nrow(x$obs)==0) return(NULL)
  cbind(site_id = s, x$obs)
}))

# adjust p-values (BH) on model-based pairwise p's
if (nrow(diffs)) {
  diffs <- diffs %>% group_by(contrast) %>% mutate(p_diff_adj = p.adjust(p_diff, method = "BH")) %>% ungroup()
}

# coords
site_parts <- unique(dat[, .(site_id, chr, pos, strand)])

# =========================
# Build wide tables for observed & predicted ratios
# =========================
ensure_cols <- function(df, cols) {
  for (cl in cols) if (!cl %in% names(df)) df[[cl]] <- NA_real_
  df
}

obs_wide <- obs_summ %>%
  dplyr::select(site_id, condition = conditions, obs_frac) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = condition, values_from = obs_frac) %>%
  ensure_cols(ref_levels) %>%
  dplyr::rename_with(~ paste0(.x, "_obs"), dplyr::all_of(ref_levels))

pred_wide <- preds %>%
  dplyr::select(site_id, condition, p_hat) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = condition, values_from = p_hat) %>%
  ensure_cols(ref_levels) %>%
  dplyr::rename_with(~ paste0(.x, "_pred"), dplyr::all_of(ref_levels))

# =========================
# Merge ratios into diffs and compute signed diffs that MATCH label A_vs_B (A - B)
# Keep both ratio sets in output
# =========================
`%||%` <- function(a,b) if (is.null(a) || is.na(a)) b else a

diffs <- diffs %>%
  tidyr::separate(contrast, into = c("A","B"), sep = "_vs_", remove = FALSE) %>%
  dplyr::left_join(site_parts, by = "site_id") %>%
  dplyr::left_join(obs_wide,  by = "site_id") %>%
  dplyr::left_join(pred_wide, by = "site_id") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    obs_A  = get(paste0(A, "_obs")),  obs_B  = get(paste0(B, "_obs")),
    pred_A = get(paste0(A, "_pred")), pred_B = get(paste0(B, "_pred")),
    obs_diff  = obs_A  - obs_B,
    pred_diff = pred_A - pred_B,
    abs_obs_diff  = abs(obs_diff),
    abs_pred_diff = abs(pred_diff)
  ) %>%
  dplyr::ungroup()

# =========================
# Sanity: observed pooling matches recomputation from long
# =========================
check_obs <- long %>%
  group_by(site_id, conditions) %>%
  summarise(ed2 = sum(ed), tot2 = sum(total), .groups = "drop") %>%
  mutate(frac2 = ifelse(tot2 > 0, ed2 / tot2, NA_real_)) %>%
  pivot_wider(names_from = conditions, values_from = frac2) %>%
  rename_with(~ paste0(.x, "_obs2"), all_of(ref_levels))

chk <- left_join(obs_wide, check_obs, by = "site_id")
if (nrow(chk)) {
  for (lv in ref_levels) {
    a <- chk[[paste0(lv, "_obs")]]
    b <- chk[[paste0(lv, "_obs2")]]
    stopifnot(all(abs(a - b) < 1e-12, na.rm = TRUE))
  }
}

# =========================
# Write outputs
# =========================

# Validate results
stopifnot(all(is.finite(diffs$obs_diff[!is.na(diffs$obs_diff)])))
stopifnot(all(is.finite(diffs$pred_diff[!is.na(diffs$pred_diff)])))
stopifnot(all(diffs$p_diff_adj >= 0 & diffs$p_diff_adj <= 1, na.rm = TRUE))

fwrite(obs_summ %>% left_join(site_parts, by="site_id"),
       sprintf("%s_observedFractions.tsv", out_prefix), sep = "\t")
fwrite(preds %>% left_join(site_parts, by="site_id"),
       sprintf("%s_predictedRatios.tsv",   out_prefix), sep = "\t")
fwrite(diffs, sprintf("%s_propDiffs.tsv", out_prefix), sep = "\t")
message("Wrote: ",
        paste(c(sprintf("%s_observedFractions.tsv", out_prefix),
                sprintf("%s_predictedRatios.tsv",   out_prefix),
                sprintf("%s_propDiffs.tsv",         out_prefix)), collapse = ", "))

# =========================
# Significant sites (FDR < 0.1 AND |obs_diff| > 0.05) and summary
# =========================
sig_sites <- diffs %>%
  dplyr::filter(is.finite(p_diff_adj), is.finite(obs_diff)) %>%
  dplyr::mutate(
    sig_FDR = p_diff_adj < FDR_THRESH,
    sig_obs = abs_obs_diff > OBS_DIFF_MIN,
    is_sig  = sig_FDR & sig_obs
  ) %>%
  dplyr::filter(is_sig)

fwrite(sig_sites, sprintf("%s_sigSites.tsv", out_prefix), sep = "\t")

all_contrasts_present <- sort(unique(diffs$contrast))
if (nrow(sig_sites)) {
  sig_summary <- sig_sites %>%
    mutate(direction = case_when(
      obs_diff > 0 ~ "A_greater_B",
      obs_diff < 0 ~ "A_less_B",
      TRUE         ~ "tie"
    )) %>%
    summarise(
      n_sig = n(),
      n_up  = sum(direction == "A_greater_B"),
      n_down= sum(direction == "A_less_B"),
      .by   = contrast
    ) %>%
    right_join(tibble(contrast = all_contrasts_present), by = "contrast") %>%
    mutate(
      n_sig = coalesce(n_sig,  0L),
      n_up  = coalesce(n_up,   0L),
      n_down= coalesce(n_down, 0L)
    ) %>%
    arrange(contrast)
} else {
  sig_summary <- tibble(
    contrast = all_contrasts_present,
    n_sig = 0L, n_up = 0L, n_down = 0L
  )
}
fwrite(sig_summary, sprintf("%s_sigSummary.tsv", out_prefix), sep = "\t")
message("Significant-site summary written: ", sprintf("%s_sigSummary.tsv", out_prefix))

# =========================
# Volcano plot: signed obs_diff vs -log10(FDR)
# not sig = grey; A>B red; A<B blue; labels use real groups
# =========================
volcano_plot_signed <- function(diffs_tbl, contrast_label,
                                num_labels = 20, outbase = "volcano_DIFF",
                                cap_y = VOLC_Y_CAP, fdr_thresh = FDR_THRESH) {
  parts <- strsplit(contrast_label, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2) { message("Bad contrast label: ", contrast_label); return(invisible(NULL)) }
  A <- parts[1]; B <- parts[2]

  dd <- diffs_tbl %>% dplyr::filter(contrast == contrast_label)
  if (!nrow(dd)) { message("No rows for contrast: ", contrast_label); return(invisible(NULL)) }

  dd <- dd %>%
    dplyr::mutate(
      x_signed  = obs_diff,  # signed (A − B)
      p_use     = dplyr::if_else(is.finite(p_diff_adj) & p_diff_adj > 0, p_diff_adj, 1),
      nlog10p_r = -log10(p_use),
      nlog10p   = pmin(nlog10p_r, cap_y),
      sig       = (p_use < fdr_thresh) & (abs(x_signed) > OBS_DIFF_MIN),
      dir_label = dplyr::case_when(
        sig & (x_signed > 0) ~ paste0(A, ">", B),
        sig & (x_signed < 0) ~ paste0(A, "<", B),
        TRUE                 ~ "not sig"
      )
    )

  desired_order  <- c("not sig", paste0(A, ">", B), paste0(A, "<", B))
  present_levels <- desired_order[desired_order %in% unique(dd$dir_label)]
  dd$dir_label   <- factor(dd$dir_label, levels = present_levels)

  labels_all <- c("not sig", paste0(A, ">", B), paste0(A, "<", B))
  colors_all <- c("grey70",  "#d62728",           "#1f77b4")
  col_map    <- setNames(colors_all[match(present_levels, labels_all)], present_levels)

  cutoff_y     <- -log10(fdr_thresh)
  cutoff_y_cap <- pmin(cutoff_y, cap_y)

  p <- ggplot2::ggplot(dd, ggplot2::aes(x = x_signed, y = nlog10p, color = dir_label)) +
    ggplot2::geom_point(alpha = 0.85, size = 1.8) +
    ggplot2::scale_color_manual(values = col_map, name = "Direction", drop = TRUE) +
    ggplot2::geom_hline(yintercept = cutoff_y_cap, linetype = 2, color = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype = 3, color = "grey75") +
    ggplot2::labs(
      title = paste0("Volcano: ", contrast_label),
      x     = paste0("observed editing ratio difference (", A, " - ", B, ")"),
      y     = expression(paste("-log"[10], "(FDR)"))
    ) +
    ggplot2::theme_minimal(base_size = 12)

  top <- dd %>% dplyr::arrange(!sig, p_use, dplyr::desc(abs(x_signed))) %>% dplyr::slice_head(n = num_labels)
  if (nrow(top)) {
    p <- p + ggrepel::geom_text_repel(
      data = top, ggplot2::aes(label = site_id),
      size = 3, max.overlaps = 50, seed = 1
    )
  }

  pdf_file <- paste0(outbase, "_", contrast_label, ".pdf")
  png_file <- paste0(outbase, "_", contrast_label, ".png")
  ggplot2::ggsave(pdf_file, p, width = 7.2, height = 6.2)
  ggplot2::ggsave(png_file, p, width = 7.2, height = 6.2, dpi = 300)
  message("Saved volcano: ", pdf_file, " and ", png_file)
}

# =========================
# Concordance plots: pred_diff vs obs_diff (signed), blue expected vs dotted fitted
# =========================
concordance_plot <- function(diffs_tbl, contrast_label,
                             outbase = "concordance") {
  parts <- strsplit(contrast_label, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2) { message("Bad contrast label: ", contrast_label); return(invisible(NULL)) }
  A <- parts[1]; B <- parts[2]

  dd <- diffs_tbl %>% dplyr::filter(contrast == contrast_label) %>%
    dplyr::select(site_id, pred_diff, obs_diff) %>%
    dplyr::filter(is.finite(pred_diff), is.finite(obs_diff))
  if (!nrow(dd)) { message("No rows for concordance: ", contrast_label); return(invisible(NULL)) }

  fit_lm    <- suppressWarnings(lm(obs_diff ~ pred_diff, data = dd))
  slope     <- unname(coef(fit_lm)[2])
  r_pearson <- suppressWarnings(cor(dd$pred_diff, dd$obs_diff, method = "pearson"))
  r_spear   <- suppressWarnings(cor(dd$pred_diff, dd$obs_diff, method = "spearman"))

  p <- ggplot(dd, aes(pred_diff, obs_diff)) +
    geom_point(alpha = 0.6, size = 1.6, color = "grey35") +
    geom_abline(aes(color = "expected", linetype = "expected"),
                slope = 1, intercept = 0, linewidth = 0.8) +
    geom_smooth(aes(color = "fitted", linetype = "fitted"),
                method = "lm", se = FALSE, linewidth = 0.9) +
    scale_color_manual(name = NULL, values = c(expected = "blue", fitted = "black"),
                       breaks = c("expected","fitted"), labels = c("expected","fitted")) +
    scale_linetype_manual(name = NULL, values = c(expected = "solid", fitted = "dotted"),
                          breaks = c("expected","fitted"), labels = c("expected","fitted")) +
    labs(
      title = paste0("Concordance: ", contrast_label),
      x = paste0("Predicted difference (", A, " − ", B, ", p̂)"),
      y = paste0("Observed difference (", A, " − ", B, ", pooled)")
    ) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.8,
             label = sprintf("slope (fitted) = %.3f\nPearson r = %.3f; Spearman \u03C1 = %.3f; n = %d",
                             slope, r_pearson, r_spear, nrow(dd))) +
    theme_minimal(base_size = 12)

  pdf_file <- paste0(outbase, "_", contrast_label, ".pdf")
  png_file <- paste0(outbase, "_", contrast_label, ".png")
  ggsave(pdf_file, p, width = 7.2, height = 6.2)
  ggsave(png_file, p, width = 7.2, height = 6.2, dpi = 300)
  message("Saved concordance: ", pdf_file, " and ", png_file)
}

# ---------- Robust list of expected one-direction contrasts (higher index minus lower index)
expected_contrasts <- apply(combn(ref_levels, 2), 2, function(p) paste0(p[2], "_vs_", p[1]))
message("Expected contrasts (one-direction): ", paste(expected_contrasts, collapse = ", "))

volc_base <- paste0(out_prefix, "_volcano_DIFF_cap", VOLC_Y_CAP)
conc_base <- paste0(out_prefix, "_concordance")

if (!nrow(diffs)) {
  message("WARNING: 'diffs' has 0 rows. No plots can be generated. Check upstream modeling/join steps.")
}

for (cc in expected_contrasts) {
  n_cc <- sum(diffs$contrast == cc, na.rm = TRUE)
  message(sprintf("Plotting %s: %d rows", cc, n_cc))
  if (n_cc > 0) {
    volcano_plot_signed(diffs, cc, outbase = volc_base, cap_y = VOLC_Y_CAP, fdr_thresh = FDR_THRESH)
    concordance_plot(diffs, cc, outbase = conc_base)
  } else {
    message("  (no rows for this contrast)")
  }
}

# =========================
# Average editing level per condition (all sites vs significant sites) + bar plot
# =========================
compute_mean_obs <- function(obs_long, site_ids = NULL, set_name = "all_sites") {
  df <- obs_long
  if (!is.null(site_ids)) df <- df %>% dplyr::filter(site_id %in% site_ids)
  df %>%
    dplyr::group_by(conditions) %>%
    dplyr::summarise(
      mean_obs = mean(obs_frac, na.rm = TRUE),
      n_sites  = dplyr::n_distinct(site_id),
      .groups = "drop"
    ) %>%
    dplyr::mutate(set = set_name)
}

obs_long <- obs_summ %>%
  dplyr::select(site_id, conditions, obs_frac) %>%
  dplyr::distinct()

sig_site_ids <- unique(sig_sites$site_id)
has_sig <- length(sig_site_ids) > 0

mean_all <- compute_mean_obs(obs_long, site_ids = NULL, set_name = "all_sites")
mean_sig <- if (has_sig) compute_mean_obs(obs_long, site_ids = sig_site_ids, set_name = "significant_sites") else {
  mean_all %>% dplyr::mutate(set = "significant_sites", mean_obs = NA_real_, n_sites = 0L)
}

mean_levels_tbl <- dplyr::bind_rows(mean_all, mean_sig) %>%
  dplyr::mutate(conditions = factor(conditions, levels = ref_levels),
                set        = factor(set, levels = c("all_sites","significant_sites")))

mean_levels_wide <- mean_levels_tbl %>%
  dplyr::select(conditions, set, mean_obs) %>%
  tidyr::pivot_wider(names_from = set, values_from = mean_obs)

data.table::fwrite(mean_levels_tbl,  sprintf("%s_meanEditing_byCondition_long.tsv", out_prefix), sep = "\t")
data.table::fwrite(mean_levels_wide, sprintf("%s_meanEditing_byCondition_wide.tsv", out_prefix), sep = "\t")

p_levels <- ggplot2::ggplot(mean_levels_tbl,
                            ggplot2::aes(x = conditions, y = mean_obs, fill = set)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_fill_manual(
    name   = "Site set",
    values = c(all_sites = "#9e9e9e", significant_sites = "#4e79a7"),
    labels = c("All sites", "Significant sites")
  ) +
  ggplot2::labs(
    title = "Average observed editing level by condition",
    x = "Condition",
    y = "Mean editing ratio"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(sprintf("%s_meanEditing_byCondition_bars.pdf", out_prefix), p_levels, width = 7.2, height = 5.0)
ggplot2::ggsave(sprintf("%s_meanEditing_byCondition_bars.png", out_prefix), p_levels, width = 7.2, height = 5.0, dpi = 300)
message("Saved mean-by-condition bar plot and tables.")

message("All done.")

# =========================
# HOW PREDICTED RATIOS & DIFFS ARE COMPUTED (reference)
# -----------------------------------------------------
# Model per site: glmmTMB with
#   cbind(ed, uned) ~ conditions (+ 1|indID if repeated)
# using beta-binomial; binomial fallback if needed.
#
# Fixed-effects linear predictor (HL baseline):
#   η_L = β0 + β_L (if L != HL), else η_HL = β0
# Predicted editing ratio per group:
#   p̂(L) = logistic(η_L)
#
# Observed editing ratio per group (pooled):
#   obs(L) = (Σ edited in L) / (Σ total in L)
#
# Single-direction contrasts only (higher index minus lower index in ref_levels):
#   LABEL = A_vs_B
#   pred_diff = p̂(A) − p̂(B)
#   obs_diff  = obs(A) − obs(B)
#   Positive values mean the labeled numerator group (A) has a larger ratio.
# =========================

