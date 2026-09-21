\name{extract_bb_contrasts}
\alias{extract_bb_contrasts}
\title{Extract coefficient and delta-method contrasts}
\description{Extracts pairwise condition contrasts from a glmmTMB fit on the coefficient scale, response-probability scale, or both.}
\usage{extract_bb_contrasts(fit, contrast_pairs, method = c("both", "coefficient", "delta"), conf_level = 0.95, engine = NULL)}
\arguments{\item{fit}{A glmmTMB model fit.}\item{contrast_pairs}{Two-column data.frame or matrix containing A and B condition labels.}\item{method}{Contrast method.}\item{conf_level}{Confidence level.}\item{engine}{Optional model engine label.}}
\value{A tibble containing the requested contrast estimates, standard errors, Wald statistics, p-values, and confidence intervals. Delta-method probability differences use Wald intervals, which are not constrained to [-1, 1].}
