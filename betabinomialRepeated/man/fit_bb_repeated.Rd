\name{fit_bb_repeated}
\alias{fit_bb_repeated}
\title{Fit beta-binomial regression with repeated measures and binomial fallback}

\usage{
fit_bb_repeated(data, site = NULL, edited = "ed", total = "total",
  condition = "conditions", subject = NULL, ref_level = NULL,
  method = c("both", "coefficient", "delta"),
  family_order = c("beta-binomial", "binomial"), min_total = 1L,
  min_samples = 2L, min_subjects = NULL, conf_level = 0.95,
  contrast_pairs = NULL, p_adjust_method = "BH",
  control = glmmTMB::glmmTMBControl(), quiet = FALSE)

}

\arguments{
\item{data}{A data.frame/tibble with one row per site-sample observation.}
\item{site}{Column name identifying sites. If NULL, all rows are treated as one site.}
\item{edited}{Column name containing edited-read counts. Must be numeric, finite, and integer-valued.}
\item{total}{Column name containing total-read counts. Must be numeric, finite, and integer-valued.}
\item{condition}{Column name containing the experimental condition.}
\item{subject}{Optional column name identifying repeated-measure subjects.}
\item{ref_level}{Reference condition. If NULL, the first observed condition is used as the reference level.}
\item{method}{One of \code{"coefficient"}, \code{"delta"}, or \code{"both"}.}
\item{family_order}{Families attempted in the specified order. The first family is the primary method; subsequent families are fallback methods if an earlier family is not accepted.}
\item{min_total}{Minimum positive total reads required for a row to enter a model. Must be at least 1.}
\item{min_samples}{Minimum number of qualifying observations across all conditions combined at a site. Filtering is not condition-specific, and a site may be retained even if it is observed in only one condition.}
\item{min_subjects}{Optional minimum number of unique subjects required across all conditions combined at a site. Ignored when \code{subject = NULL}.}
\item{conf_level}{Confidence level for intervals.}
\item{contrast_pairs}{Optional two-column data.frame/matrix defining comparisons. By default all unordered pairs are generated and reported as \code{A_vs_B} with the first group minus the second group. In site-wise analyses, a requested level that is unavailable at a particular site produces a structured site-level \code{not_contrastable} result rather than aborting the entire analysis.}
\item{p_adjust_method}{Multiple-testing adjustment method passed to \code{stats::p.adjust()}. Default is \code{"BH"}. Use \code{"none"} to retain raw p-values. Adjustment is performed separately within each pairwise contrast and separately for coefficient- and delta-method p-values.}
\item{control}{A \code{glmmTMBControl()} object controlling model optimization.}
\item{quiet}{Suppress per-site messages.}
}
\description{
Fits one model per site using a beta-binomial mixed model when possible, with automatic binomial fallback. Repeated measures are represented by a random intercept for \code{subject}.
\cr
\code{method = "coefficient"} returns contrasts on the model's log-odds scale. \code{method = "delta"} returns contrasts of fitted probabilities using the delta method. When a subject random intercept is present, these are fixed-effect probabilities evaluated at a random effect of zero, not population-averaged (marginal) probabilities.
\cr
By default, raw Wald p-values are adjusted using the Benjamini-Hochberg (BH) procedure separately within each pairwise contrast and separately for coefficient- and delta-method p-values.
}
\value{An object of class \code{bb_repeated}. Each element of \code{fits} has a \code{status} of \code{success}, \code{model_failed}, or \code{not_contrastable}. The top-level object reports \code{n_sites_retained}, \code{n_sites_attempted}, \code{n_sites_successful}, \code{n_sites_failed}, and \code{n_sites_not_contrastable}.}
\examples{
\dontrun{
fit <- fit_bb_repeated(dat, site="site_id", edited="ed", total="total",
                       condition="conditions", subject="indID",
                       ref_level="HL", method="both")
fit$fits[["chr1:123:+"]]$contrasts
fit$fits[["chr1:123:+"]]$comparison
}
}

