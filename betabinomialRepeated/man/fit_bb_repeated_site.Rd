\name{fit_bb_repeated_site}
\alias{fit_bb_repeated_site}
\title{Fit one site-wise beta-binomial or binomial model}
\description{Fits one site and returns the selected model, diagnostics, and contrast results.}
\usage{fit_bb_repeated_site(data, condition_col = ".condition__", edited_col = ".edited__", total_col = ".total__", subject_col = NULL, method = c("both", "coefficient", "delta"), family_order = c("beta-binomial", "binomial"), conf_level = 0.95, contrast_pairs = NULL, control = glmmTMB::glmmTMBControl(), quiet = TRUE)}
\arguments{\item{data}{Prepared site-level data.}\item{condition_col}{Condition column name.}\item{edited_col}{Edited count column name.}\item{total_col}{Total count column name.}\item{subject_col}{Optional subject column name.}\item{method}{Contrast scale to return.}\item{family_order}{Families attempted in order.}\item{conf_level}{Confidence level.}\item{contrast_pairs}{Optional pairwise contrasts. If a requested level is unavailable at the site, a structured failure result is returned before model fitting.}\item{control}{glmmTMB optimizer control object.}\item{quiet}{Suppress messages.}}
\value{A list containing the fitted model, engine, contrasts, comparison table, and fitting errors.}
