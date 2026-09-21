\name{compare_bb_methods}
\alias{compare_bb_methods}
\title{Compare coefficient and delta-method contrasts}
\description{Compares coefficient-scale and response-scale contrast estimates and their statistical conclusions. Raw and adjusted significance indicators are reported separately. When adjusted p-values are available, \code{significance_agreement} uses the adjusted indicators; otherwise it uses the raw indicators.}
\usage{compare_bb_methods(contrasts)}
\arguments{\item{contrasts}{Output from \code{extract_bb_contrasts} with both methods.}}
\value{A tibble containing absolute estimate differences, direction agreement, raw significance indicators, adjusted significance indicators when available, and primary significance agreement.}
