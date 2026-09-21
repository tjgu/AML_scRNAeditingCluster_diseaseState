#!/usr/bin/env Rscript
# File: association_editingGene_betaBinomial_nestModel.R
# Purpose:
#   Beta-binomial regression of global RNA editing (sum_ed/sum_total)
#   against each gene's expression, adjusting for condition (HL/ND/RM/PO).
#   Also fits interaction (expr * condition) and provides:
#     - overall association adjusted for condition
#     - interaction test (does slope differ by condition?)
#     - interaction term coefficients per non-baseline condition (which condition differs?)
#     - per-condition slope (estimated marginal trend) Wald tests with BH-FDR
#
# Outputs:
#   - overall_betabinom_gene_assoc.tsv.gz
#   - overall_betabinom_gene_assoc_FDRsig.tsv
#   - interaction_betabinom_gene_assoc.tsv.gz
#   - interaction_betabinom_gene_assoc_FDRsig.tsv
#   - perCondition_slopes_from_interaction.tsv.gz
#   - perCondition_wald_slopes.tsv.gz
#   - perCondition_wald_slopes_FDRsig.tsv
#   - summary.txt
#
# Notes:
#   Requires: data.table, glmmTMB
#   Models:
#     M0: cbind(sum_ed, sum_total - sum_ed) ~ conditions
#     M1: cbind(sum_ed, sum_total - sum_ed) ~ expr + conditions
#     M2: cbind(sum_ed, sum_total - sum_ed) ~ expr * conditions
#   expr = log2(FPKM + 0.1) by default
#
# Changes requested:
#   - Gene prevalence filter uses expression > 3
#   - Interaction table includes per-condition interaction terms (ND/RM/PO vs baseline HL) with BH-FDR

suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
})

options(stringsAsFactors = FALSE)

# =========================
# CONFIG
# =========================
GENE_EXPR_FILE <- "./data/rsem_fpkm_matrix.withGene.tsv"
META_FILE      <- "./data/meta_samples.txt"
EDIT_FILE      <- "./results/de_bulk_sites/de_bulk_sites_filteredLong.tsv.gz"

OUT_DIR <- paste0("./results/editingAllGenes_betaBinom_", format(Sys.Date(), "%Y%m%d"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

alpha_fdr       <- 0.1
min_expr_prop   <- 0.25
expr_threshold  <- 3
log_expr        <- TRUE
min_n_samples   <- 10
conditions_level_order <- c("HL","ND","RM","PO")

do_interaction <- TRUE  # required for per-condition tests

# =========================
# Helpers
# =========================
normalize_seqid <- function(x) {
  x <- as.character(x)
  sub("^X(?=\\d)", "", x, perl = TRUE)
}
num <- function(x) suppressWarnings(as.numeric(x))

fit_bb <- function(df, formula) {
  tryCatch(
    withCallingHandlers(
      glmmTMB(formula, data = df, family = betabinomial(link = "logit")),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
}

get_coef_row <- function(m, term) {
  if (is.null(m)) return(NULL)
  sm <- summary(m)
  coefs <- sm$coefficients$cond
  if (is.null(coefs) || !(term %in% rownames(coefs))) return(NULL)
  out <- as.list(coefs[term, , drop = TRUE])
  data.table(
    estimate  = unname(out[["Estimate"]]),
    std_error = unname(out[["Std. Error"]]),
    statistic = unname(out[["z value"]]),
    p_value   = unname(out[["Pr(>|z|)"]])
  )
}

lrt_p <- function(m0, m1) {
  if (is.null(m0) || is.null(m1)) return(NA_real_)
  aa <- tryCatch(anova(m0, m1), error = function(e) NULL)
  if (is.null(aa) || nrow(aa) < 2) return(NA_real_)
  pcol <- intersect(c("Pr(>Chisq)"), colnames(aa))
  if (!length(pcol)) return(NA_real_)
  as.numeric(aa[2, pcol[1], drop = TRUE])
}

# Extract interaction term stats for each non-baseline condition from interaction model m2
# Returns a named list of columns: int_est_<cond>, int_se_<cond>, int_z_<cond>, int_p_<cond>
get_interaction_terms <- function(m2, cond_lvls, baseline) {
  out <- list()
  for (cc in cond_lvls) {
    if (cc == baseline) next
    out[[paste0("int_est_", cc)]] <- NA_real_
    out[[paste0("int_se_",  cc)]] <- NA_real_
    out[[paste0("int_z_",   cc)]] <- NA_real_
    out[[paste0("int_p_",   cc)]] <- NA_real_
  }
  if (is.null(m2)) return(out)

  sm <- summary(m2)
  coefs <- sm$coefficients$cond
  if (is.null(coefs)) return(out)

  find_term <- function(cc) {
    t1 <- paste0("expr:conditions", cc)
    t2 <- paste0("conditions", cc, ":expr")
    if (t1 %in% rownames(coefs)) return(t1)
    if (t2 %in% rownames(coefs)) return(t2)
    return(NA_character_)
  }

  for (cc in cond_lvls) {
    if (cc == baseline) next
    term <- find_term(cc)
    if (is.na(term)) next
    out[[paste0("int_est_", cc)]] <- coefs[term, "Estimate"]
    out[[paste0("int_se_",  cc)]] <- coefs[term, "Std. Error"]
    out[[paste0("int_z_",   cc)]] <- coefs[term, "z value"]
    out[[paste0("int_p_",   cc)]] <- coefs[term, "Pr(>|z|)"]
  }
  out
}

# Per-condition slope + Wald test from interaction model:
# slope(cond) = b_expr + b_expr:cond (baseline has only b_expr)
# var(slope)  = var(b_expr) + var(b_int) + 2*cov(b_expr, b_int)
extract_condition_slopes_wald <- function(m2, cond_lvls, baseline) {
  if (is.null(m2)) return(NULL)

  b <- tryCatch(fixef(m2)$cond, error = function(e) NULL)
  V <- tryCatch(vcov(m2)$cond,  error = function(e) NULL)
  if (is.null(b) || is.null(V) || is.null(names(b))) return(NULL)

  int_term_name <- function(cc) {
    t1 <- paste0("expr:conditions", cc)
    t2 <- paste0("conditions", cc, ":expr")
    if (t1 %in% names(b)) return(t1)
    if (t2 %in% names(b)) return(t2)
    return(NA_character_)
  }

  if (!("expr" %in% names(b))) return(NULL)
  b_expr <- b[["expr"]]
  v_expr <- V["expr","expr"]

  out <- vector("list", length(cond_lvls))

  for (k in seq_along(cond_lvls)) {
    cc <- cond_lvls[k]
    if (cc == baseline) {
      slope <- b_expr
      se <- sqrt(v_expr)
    } else {
      it <- int_term_name(cc)
      if (is.na(it)) {
        slope <- NA_real_
        se <- NA_real_
      } else {
        slope <- b_expr + b[[it]]
        v_it <- V[it, it]
        c_xy <- V["expr", it]
        se <- sqrt(v_expr + v_it + 2*c_xy)
      }
    }

    z <- if (is.finite(slope) && is.finite(se) && se > 0) slope / se else NA_real_
    p <- if (is.finite(z)) 2 * pnorm(-abs(z)) else NA_real_

    out[[k]] <- data.table(
      conditions = cc,
      slope_logit = slope,
      slope_se = se,
      z = z,
      p_value_wald = p
    )
  }

  rbindlist(out)
}

# =========================
# Load metadata
# =========================
meta <- fread(META_FILE, na.strings = c("NA","NaN","nan",""))
low <- tolower(names(meta))
col_seq <- if ("seqid" %in% low) names(meta)[which(low=="seqid")[1]] else if ("sample" %in% low) names(meta)[which(low=="sample")[1]] else stop("seqID column not found in meta")
col_con <- if ("conditions" %in% low) names(meta)[which(low=="conditions")[1]] else if ("condition" %in% low) names(meta)[which(low=="condition")[1]] else stop("conditions column not found in meta")
col_ind <- if ("indid" %in% low) names(meta)[which(low=="indid")[1]] else NA_character_

setnames(meta, old = col_seq, new = "seqID")
setnames(meta, old = col_con, new = "conditions")
if (!is.na(col_ind)) setnames(meta, old = col_ind, new = "indID") else meta[, indID := NA_character_]

meta[, `:=`(
  seqID = normalize_seqid(seqID),
  conditions = as.character(conditions),
  indID = as.character(indID)
)]

# =========================
# Load editing long; aggregate global counts per sample
# =========================
edit_long <- fread(EDIT_FILE, na.strings = c("NA","NaN","nan",""))
if (ncol(edit_long) < 7) stop("Editing long file must have at least 7 columns per provided spec.")

idx <- unique(c(1, 2, 3, 4, 5, 6, 7, ncol(edit_long)))
std_names <- c("edit_id","seqID","uned","ed","total","present","conditions","indID")
setnames(edit_long, old = names(edit_long)[idx], new = std_names)

edit_long[, `:=`(
  seqID = normalize_seqid(seqID),
  conditions = as.character(conditions),
  ed = num(ed), uned = num(uned), total = num(total),
  present = as.integer(present), indID = as.character(indID)
)]

agg_by_sample <- edit_long[is.finite(ed) & is.finite(total) & total > 0,
  .(sum_ed = sum(ed, na.rm = TRUE),
    sum_total = sum(total, na.rm = TRUE),
    ratio_weighted = sum(ed, na.rm = TRUE) / sum(total, na.rm = TRUE),
    n_sites = .N),
  by = .(seqID)
]

per_sample <- merge(agg_by_sample, unique(meta[, .(seqID, conditions, indID)]),
                    by = "seqID", all.x = TRUE)

extra_lvls <- setdiff(sort(unique(as.character(per_sample$conditions))), conditions_level_order)
per_sample[, conditions := factor(conditions, levels = c(conditions_level_order, extra_lvls))]

per_sample <- per_sample[!is.na(conditions) & is.finite(sum_ed) & is.finite(sum_total) & sum_total > 0]
per_sample <- per_sample[sum_ed >= 0 & sum_ed <= sum_total]

# =========================
# Load gene expression matrix
# =========================
expr_dt <- fread(GENE_EXPR_FILE, na.strings = c("NA","NaN","nan",""))
if (ncol(expr_dt) < 3) stop("Gene expression matrix must have at least 3 columns (gene_id, gene_name, samples)")
setnames(expr_dt, old = names(expr_dt)[1:2], new = c("gene_id","gene_name"))
expr_dt[, gene_name := as.character(gene_name)]

sample_cols_raw  <- setdiff(names(expr_dt), c("gene_id","gene_name"))
sample_cols_norm <- normalize_seqid(sample_cols_raw)

if (any(duplicated(sample_cols_norm))) {
  dup <- unique(sample_cols_norm[duplicated(sample_cols_norm)])
  stop(sprintf("Duplicate sample IDs after normalization: %s", paste(dup, collapse = ", ")))
}
setnames(expr_dt, old = sample_cols_raw, new = sample_cols_norm)
sample_cols <- sample_cols_norm

# =========================
# Harmonize samples
# =========================
samples_use <- intersect(sample_cols, per_sample$seqID)
if (length(samples_use) < min_n_samples) stop("Not enough overlapping samples between expression and editing data.")

ps_use <- per_sample[match(samples_use, per_sample$seqID)]
stopifnot(all(ps_use$seqID == samples_use))

expr_sub <- expr_dt[, c("gene_id","gene_name", samples_use), with = FALSE]
expr_sub[, (samples_use) := lapply(.SD, num), .SDcols = samples_use]

# =========================
# Filter genes by prevalence
# =========================
n_samples <- length(samples_use)
min_nonzero <- max(3L, ceiling(min_expr_prop * n_samples))
nz_counts <- expr_sub[, rowSums(as.matrix(.SD) > expr_threshold, na.rm = TRUE), .SDcols = samples_use]
keep_idx  <- nz_counts >= min_nonzero
expr_sub  <- expr_sub[keep_idx]

message(sprintf(
  "Kept %d genes after prevalence filter: expression > %s in >= %.1f%% samples (>= %d of %d).",
  nrow(expr_sub), as.character(expr_threshold), 100*min_expr_prop, min_nonzero, n_samples
))

E <- as.matrix(expr_sub[, ..samples_use])
if (log_expr) E <- log2(E + 0.1)

base_df <- data.table(
  seqID      = ps_use$seqID,
  conditions = ps_use$conditions,
  sum_ed     = as.integer(round(ps_use$sum_ed)),
  sum_total  = as.integer(round(ps_use$sum_total))
)
base_df[, sum_uned := sum_total - sum_ed]
base_df <- base_df[sum_uned >= 0]

# =========================
# Per-gene beta-binomial regression
# =========================
run_one_gene <- function(i) {
  gene_id   <- expr_sub$gene_id[i]
  gene_name <- expr_sub$gene_name[i]
  x <- as.numeric(E[i, ])

  df <- copy(base_df)
  df[, expr := x]
  df <- df[is.finite(expr) & !is.na(conditions)]
  if (nrow(df) < min_n_samples) return(NULL)

  # M0: no expr
  m0 <- fit_bb(df, cbind(sum_ed, sum_uned) ~ conditions)
  # M1: expr + conditions
  m1 <- fit_bb(df, cbind(sum_ed, sum_uned) ~ expr + conditions)

  coef_expr <- get_coef_row(m1, "expr")
  if (is.null(coef_expr)) {
    coef_expr <- data.table(estimate=NA_real_, std_error=NA_real_, statistic=NA_real_, p_value=NA_real_)
  }

  p_lrt <- lrt_p(m0, m1)

  out <- data.table(
    gene_id = gene_id,
    gene_name = gene_name,
    n_used = nrow(df),
    estimate = coef_expr$estimate,
    std_error = coef_expr$std_error,
    statistic = coef_expr$statistic,
    p_value_wald = coef_expr$p_value,
    p_value_lrt  = p_lrt
  )

  if (!isTRUE(do_interaction)) return(out)

  # M2: expr * conditions
  m2 <- fit_bb(df, cbind(sum_ed, sum_uned) ~ expr * conditions)
  p_int <- lrt_p(m1, m2)
  out[, p_interaction_lrt := p_int]

  # Add per-condition interaction term coefficients (ND/RM/PO vs baseline)
  cond_lvls <- levels(df$conditions)
  baseline <- cond_lvls[1]
  int_terms <- get_interaction_terms(m2, cond_lvls, baseline)
  out[, names(int_terms) := int_terms]

  list(main = out, m_int = m2)
}

message(sprintf("Fitting beta-binomial models for %d genes ...", nrow(expr_sub)))

res_main <- vector("list", nrow(expr_sub))
res_intm <- vector("list", nrow(expr_sub))

for (i in seq_len(nrow(expr_sub))) {
  rr <- run_one_gene(i)
  if (isTRUE(do_interaction)) {
    if (is.list(rr) && !is.null(rr$main)) {
      res_main[[i]] <- rr$main
      res_intm[[i]] <- rr$m_int
    }
  } else {
    res_main[[i]] <- rr
  }
  if (i %% 500 == 0) message(sprintf("  ... %d / %d", i, nrow(expr_sub)))
}

overall_dt <- rbindlist(res_main, use.names = TRUE, fill = TRUE)
overall_dt[, fdr_wald := p.adjust(p_value_wald, method = "BH")]
overall_dt[, fdr_lrt  := p.adjust(p_value_lrt,  method = "BH")]
setorder(overall_dt, fdr_lrt, fdr_wald, p_value_lrt, p_value_wald)

fwrite(overall_dt,
       file = file.path(OUT_DIR, "overall_betabinom_gene_assoc.tsv.gz"),
       sep = "\t")
fwrite(overall_dt[fdr_lrt <= alpha_fdr],
       file = file.path(OUT_DIR, "overall_betabinom_gene_assoc_FDRsig.tsv"),
       sep = "\t")

# =========================
# Interaction outputs + per-condition slopes
# =========================
int_dt <- NULL
slopes_dt <- NULL
percond_wald_dt <- NULL

if (isTRUE(do_interaction)) {
  # Build interaction table including omnibus interaction p/FDR + interaction terms
  cond_lvls <- levels(base_df$conditions)
  baseline <- cond_lvls[1]
  nonbase <- cond_lvls[cond_lvls != baseline]

  int_cols <- c("gene_id","gene_name","n_used","estimate","std_error","statistic",
                "p_value_wald","p_value_lrt","p_interaction_lrt")

  for (cc in nonbase) {
    for (suffix in c("est","se","z","p")) {
      nm <- paste0("int_", suffix, "_", cc)
      if (nm %in% names(overall_dt)) int_cols <- c(int_cols, nm)
    }
  }

  int_dt <- overall_dt[, ..int_cols]
  int_dt[, fdr_interaction := p.adjust(p_interaction_lrt, method = "BH")]

  # Add BH-FDR for each interaction term p-value across genes
  for (cc in nonbase) {
    pnm <- paste0("int_p_", cc)
    fnm <- paste0("int_fdr_", cc)
    if (pnm %in% names(int_dt)) {
      int_dt[, (fnm) := p.adjust(get(pnm), method = "BH")]
    }
  }

  setorder(int_dt, fdr_interaction, p_interaction_lrt)

  fwrite(int_dt,
         file = file.path(OUT_DIR, "interaction_betabinom_gene_assoc.tsv.gz"),
         sep = "\t")
  fwrite(int_dt[fdr_interaction <= alpha_fdr],
         file = file.path(OUT_DIR, "interaction_betabinom_gene_assoc_FDRsig.tsv"),
         sep = "\t")

  # --- Slopes only (legacy) ---
  slope_list <- list()

  # --- Per-condition Wald tests table ---
  wald_list <- list()

  for (i in seq_along(res_intm)) {
    m2 <- res_intm[[i]]
    if (is.null(m2)) next

    gene_id   <- expr_sub$gene_id[i]
    gene_name <- expr_sub$gene_name[i]

    sm <- summary(m2)
    coefs <- sm$coefficients$cond
    if (!is.null(coefs)) {
      b_expr <- if ("expr" %in% rownames(coefs)) coefs["expr","Estimate"] else NA_real_

      for (cc in cond_lvls) {
        slope <- b_expr
        if (!is.na(slope) && cc != baseline) {
          term1 <- paste0("expr:conditions", cc)
          term2 <- paste0("conditions", cc, ":expr")
          if (term1 %in% rownames(coefs)) slope <- slope + coefs[term1,"Estimate"]
          if (term2 %in% rownames(coefs)) slope <- slope + coefs[term2,"Estimate"]
        }
        slope_list[[length(slope_list) + 1]] <- data.table(
          gene_id = gene_id,
          gene_name = gene_name,
          conditions = cc,
          slope_logit = slope
        )
      }
    }

    tmp <- extract_condition_slopes_wald(m2, cond_lvls, baseline)
    if (!is.null(tmp) && nrow(tmp)) {
      tmp[, `:=`(gene_id = gene_id, gene_name = gene_name)]
      setcolorder(tmp, c("gene_id","gene_name","conditions","slope_logit","slope_se","z","p_value_wald"))
      wald_list[[length(wald_list) + 1]] <- tmp
    }
  }

  slopes_dt <- rbindlist(slope_list, use.names = TRUE, fill = TRUE)
  setorder(slopes_dt, gene_name, conditions)
  fwrite(slopes_dt,
         file = file.path(OUT_DIR, "perCondition_slopes_from_interaction.tsv.gz"),
         sep = "\t")

  percond_wald_dt <- rbindlist(wald_list, use.names = TRUE, fill = TRUE)
  if (!is.null(percond_wald_dt) && nrow(percond_wald_dt)) {
    percond_wald_dt[, fdr_wald := p.adjust(p_value_wald, method = "BH"), by = .(conditions)]
    setorder(percond_wald_dt, conditions, fdr_wald, p_value_wald, gene_name)

    fwrite(percond_wald_dt,
           file = file.path(OUT_DIR, "perCondition_wald_slopes.tsv.gz"),
           sep = "\t")

    fwrite(percond_wald_dt[fdr_wald <= alpha_fdr],
           file = file.path(OUT_DIR, "perCondition_wald_slopes_FDRsig.tsv"),
           sep = "\t")
  }
}

# =========================
# Summary (include pos/neg significant counts overall + per condition)
# =========================
sig_overall <- overall_dt[fdr_lrt <= alpha_fdr & is.finite(estimate)]
n_sig_overall <- nrow(sig_overall)
n_pos_overall <- sum(sig_overall$estimate > 0, na.rm = TRUE)
n_neg_overall <- sum(sig_overall$estimate < 0, na.rm = TRUE)

summary_lines <- c(
  sprintf("Samples used: %d", nrow(base_df)),
  sprintf("Genes tested: %d", nrow(expr_sub)),
  sprintf("Overall (add expr) hits at FDR_lrt<=%.3f: %d", alpha_fdr, n_sig_overall),
  sprintf("  Positive (estimate>0): %d", n_pos_overall),
  sprintf("  Negative (estimate<0): %d", n_neg_overall),
  sprintf("Gene prevalence filter: expression > %s in >= %.1f%% samples", as.character(expr_threshold), 100*min_expr_prop)
)

if (isTRUE(do_interaction) && !is.null(int_dt)) {
  sig_int <- int_dt[fdr_interaction <= alpha_fdr]
  n_sig_int <- nrow(sig_int)

  sig_int2 <- merge(sig_int[, .(gene_id, gene_name, fdr_interaction)],
                    overall_dt[, .(gene_id, estimate)], by = c("gene_id"), all.x = TRUE)
  sig_int2 <- sig_int2[is.finite(estimate)]
  n_pos_int <- sum(sig_int2$estimate > 0, na.rm = TRUE)
  n_neg_int <- sum(sig_int2$estimate < 0, na.rm = TRUE)

  summary_lines <- c(summary_lines,
    sprintf("Interaction (expr:condition) hits at FDR_interaction<=%.3f: %d", alpha_fdr, n_sig_int),
    sprintf("  Positive (overall estimate>0): %d", n_pos_int),
    sprintf("  Negative (overall estimate<0): %d", n_neg_int)
  )
}

if (!is.null(percond_wald_dt) && nrow(percond_wald_dt)) {
  summary_lines <- c(summary_lines, "Per-condition Wald slope tests (from interaction model):")
  for (cc in levels(base_df$conditions)) {
    dtc <- percond_wald_dt[conditions == cc & is.finite(slope_logit)]
    sigc <- dtc[fdr_wald <= alpha_fdr]
    n_sig <- nrow(sigc)
    n_pos <- sum(sigc$slope_logit > 0, na.rm = TRUE)
    n_neg <- sum(sigc$slope_logit < 0, na.rm = TRUE)
    summary_lines <- c(summary_lines,
      sprintf("  %s: sig=%d (pos=%d, neg=%d) at FDR<=%.3f", cc, n_sig, n_pos, n_neg, alpha_fdr)
    )
  }
}

writeLines(summary_lines, con = file.path(OUT_DIR, "summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
cat("\nOutputs written to:", normalizePath(OUT_DIR, mustWork = FALSE), "\n", sep = " ")
cat(" - overall_betabinom_gene_assoc.tsv.gz\n")
cat(" - overall_betabinom_gene_assoc_FDRsig.tsv\n")
if (isTRUE(do_interaction)) {
  cat(" - interaction_betabinom_gene_assoc.tsv.gz\n")
  cat(" - interaction_betabinom_gene_assoc_FDRsig.tsv\n")
  cat(" - perCondition_slopes_from_interaction.tsv.gz\n")
  cat(" - perCondition_wald_slopes.tsv.gz\n")
  cat(" - perCondition_wald_slopes_FDRsig.tsv\n")
}
cat(" - summary.txt\n")

