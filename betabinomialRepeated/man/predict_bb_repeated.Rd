\name{predict_bb_repeated}
\alias{predict_bb_repeated}
\title{Predict condition-specific editing probabilities}
\usage{
predict_bb_repeated(fit, newdata = NULL,
                     type = c("conditional", "marginal", "both"))
}
\arguments{
\item{fit}{A glmmTMB fit.}
\item{newdata}{Optional data frame for prediction. If omitted, predictions are generated for every condition level in the fitted model.}
\item{type}{Prediction scale: "conditional" (fixed-effect probability at random effect zero), "marginal" (integrated over the fitted random-intercept distribution), or "both".}
}
\value{
A tibble containing predicted probabilities.
}
\description{
The fitted regression model is unchanged. When a subject random intercept is present, it remains in the glmmTMB model used to estimate the fixed effects and their covariance. Conditional predictions use the fixed-effect linear predictor with the random effect set to zero. Marginal predictions numerically integrate the conditional response probability over the estimated normal random-intercept distribution. Marginal predictions are descriptive population-averaged probabilities; this function does not calculate marginal-effect standard errors or confidence intervals.
}
\details{
With type = "conditional", the returned column is "predicted_probability". With type = "marginal", the same column contains the marginal probability. With type = "both", the returned columns are "conditional_probability" and "marginal_probability".
}
\keyword{models}
