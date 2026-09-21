# scRNA-seq RNA Editing Analysis in Acute Myeloid Leukemia (AML)

This repository contains R scripts and the `betabinomialRepeated` R package used to analyze RNA editing from single-cell RNA-sequencing (scRNA-seq) data across AML disease states.

The analysis uses sample-level pseudo-bulk editing counts and beta-binomial models to account for overdispersion, with repeated-measures random effects where applicable.

## Analysis Workflow

The repository contains three main components.

### 1. Differential RNA Editing Analysis

**`de_bulkEditing_betaBinomial_binomial_repeatedMeasure.R`**

This script tests differences in RNA editing at individual sites across AML disease states.

Key features include:

* Sample-level pseudo-bulk editing counts.
* Beta-binomial regression for overdispersed editing counts.
* Binomial regression fallback when the beta-binomial model cannot be used.
* Individual-level random intercepts for repeated measurements.
* Pairwise comparisons between disease states.
* Delta-method inference for differences in predicted editing probabilities.
* Benjamini-Hochberg FDR correction.
* Effect-size filtering and visualization, including volcano and concordance plots.

The primary model is:

```r
cbind(edited, unedited) ~ condition + (1 | indID)
```

### 2. Gene-Level Editing–Expression Association

**`association_editingGene_betaBinomial_nestModel.R`**

This script evaluates associations between gene expression and global RNA editing across AML disease states using nested beta-binomial models.

For each gene, the analysis evaluates:

* Overall association between gene expression and global editing.
* Expression-by-disease-state interaction.
* Condition-specific expression–editing slopes.

Gene expression is modeled as:

```r
log2(FPKM + 0.1)
```

Multiple-testing correction is performed using the Benjamini-Hochberg procedure.

### 3. `betabinomialRepeated` R Package

**`betabinomialRepeated/`**

`betabinomialRepeated` is a reusable R package for beta-binomial regression with optional subject-level repeated-measure random effects.

The package provides:

* Beta-binomial regression with binomial fallback.
* Subject-level random intercepts for repeated measurements.
* Pairwise condition contrasts.
* Coefficient-scale and probability-scale inference.
* Delta-method contrasts.
* Conditional and marginal predictions.
* Site-level model diagnostics and status reporting.

The package is maintained within this repository so that the statistical framework used for the analysis can be reused and developed alongside the analysis scripts.

## Statistical Framework

RNA editing is represented using edited and unedited read counts rather than modeling editing ratios directly:

```r
cbind(edited, unedited)
```

The beta-binomial model accounts for overdispersion relative to the standard binomial model. For repeated-measure analyses, an individual-level random intercept accounts for correlation among samples from the same individual.

For differential editing, condition-specific editing probabilities and pairwise differences are estimated from the fitted models. Multiple hypothesis testing is controlled using the Benjamini-Hochberg procedure.

For the current AML site-level analysis, the condition order is:

```text
HL, ND, RM, PO
```

Pairwise contrasts are reported as **A − B**.

## Analysis Overview

```text
scRNA-seq RNA editing data
          │
          ▼
   Sample-level
   pseudo-bulk counts
          │
          ├──────────────────────────────┐
          ▼                              ▼
 Differential editing             Gene-level association
      analysis                         analysis
          │                              │
          ▼                              ▼
 Beta-binomial                    Nested beta-binomial
 repeated-measures models             models
          │                              │
          ▼                              ▼
 Pairwise editing                Expression–editing
   comparisons                     associations
          │                              │
          └──────────────┬───────────────┘
                         ▼
                  FDR-controlled
                    results
```

## Repository Structure

```text
AML_scRNAeditingCluster_diseaseState/
├── README.md
├── de_bulkEditing_betaBinomial_binomial_repeatedMeasure.R
├── association_editingGene_betaBinomial_nestModel.R
└── betabinomialRepeated/
    ├── DESCRIPTION
    ├── NAMESPACE
    ├── R/
    ├── man/
    ├── tests/
    ├── inst/
    └── README.md
```

## Requirements

The analysis scripts require R and packages including:

* `data.table`
* `dplyr`
* `glmmTMB`
* `tidyr`
* `broom.mixed`
* `ggplot2`
* `ggrepel`
* `scales`

The `betabinomialRepeated` package requires R ≥ 4.1.0 and depends on `glmmTMB`, `dplyr`, and `tibble`.

## Package Version

`betabinomialRepeated` version **0.1.0**.

## Purpose

The goal of this repository is to provide a reproducible computational framework for studying RNA editing patterns across AML disease states and their association with gene expression using scRNA-seq-derived pseudo-bulk data.
