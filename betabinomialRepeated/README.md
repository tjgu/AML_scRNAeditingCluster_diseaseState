# betabinomialRepeated

An R package for site-level RNA editing analysis using beta-binomial or binomial regression, with optional repeated-measures random effects, automatic fallback between model families, pairwise condition contrasts, delta-method inference, and optional marginal prediction.

## Overview

`betabinomialRepeated` models edited and unedited read counts for individual RNA editing sites. Beta-binomial regression is supported for overdispersed count data, with binomial regression available as an alternative or automatic fallback according to the user-specified `family_order`.

When repeated measurements are available, the model includes a subject-specific random intercept:

```r
cbind(edited, unedited) ~ condition + (1 | subject)
```

The random effect accounts for correlation among repeated observations from the same subject.

## Site filtering

`min_samples` and `min_subjects` are evaluated at the site level across all conditions combined. Sites do not need to be represented in every condition to be retained.

A site represented in only one condition is retained but classified as `not_contrastable`, because no between-condition contrast can be estimated for that site.

## Model selection and fallback

The first family in `family_order` is the primary model. Subsequent families are fallback methods.

For example:

```r
family_order = c("beta-binomial", "binomial")
```

uses beta-binomial regression as the primary method and binomial regression as the fallback.

Conversely:

```r
family_order = c("binomial", "beta-binomial")
```

uses binomial regression as the primary method.

Diagnostics are retained separately for each attempted family.

## Contrasts

Pairwise condition contrasts can be supplied explicitly. A contrast that is unavailable for one site does not interrupt the analysis of other sites. Such a site is retained with status `not_contrastable`.

## Delta-method inference

For repeated-measures models, delta-method probabilities are calculated from the fixed-effect component of the fitted model:

```r
p = plogis(X %*% beta)
```

These are conditional probabilities evaluated at a random effect of zero (`b = 0`). They are not population-averaged marginal probabilities.

## Marginal prediction

An optional marginal prediction integrates the conditional probability over the estimated random-intercept distribution. This provides a population-averaged descriptive alternative to the primary `b = 0` conditional prediction.

## Multiple testing

Raw and adjusted p-values are reported separately. By default, BH adjustment
is performed separately within each pairwise contrast and separately for
coefficient- and delta-method p-values. Thus, for example, all site-level
p-values for `A_vs_B` form one multiple-testing family, while `A_vs_C` forms
a separate family. The primary `significance_agreement` indicator uses
adjusted significance when adjusted p-values are available.

## Input validation

`edited` and `total` must be numeric, finite, nonnegative, integer-valued count variables, with:

```r
edited <= total
```

Factor and character count columns are rejected rather than silently converted.

## Site status

Each retained site is assigned one of:

- `success`: model fitted and requested contrasts are available.
- `model_failed`: model fitting was attempted but no requested model family was accepted.
- `not_contrastable`: the site cannot produce the requested between-condition contrast, including sites represented in only one condition.

Top-level results report:

- `n_sites_retained`
- `n_sites_attempted`
- `n_sites_successful`
- `n_sites_failed`
- `n_sites_not_contrastable`
