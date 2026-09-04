# Changelog

## **postlink 0.1.2**

#### Bug Fixes

- Migrated `.stan` models to the new Stan `array[...]` syntax to resolve
  compiler parsing errors on R-devel and platforms running `stanc3`
  compilers.
- Increased minimum dependency requirements for `rstan` and
  `StanHeaders` to `>= 2.26.0`.

## **postlink 0.1.1**

#### New Features

- Updated the `priors` argument for the Bayesian Mixture functions
  ([`adjMixBayes()`](https://postlink-group.github.io/postlink/reference/adjMixBayes.md),
  [`glmMixBayes()`](https://postlink-group.github.io/postlink/reference/glmMixBayes.md),
  [`survregMixBayes()`](https://postlink-group.github.io/postlink/reference/survregMixBayes.md))
  so that intercept priors and slope coefficient priors are now
  decoupled.

#### Architecture & S3 Methods

- Refactored S3 class assignment system to remove false inheritance from
  base `glm`, `lm`, and `coxph` classes. Fitted objects now use
  dedicated package-level classes (e.g., `plmodel`, `plglm`, `plcoxph`).
- Added standard generic extractors (i.e.,
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`df.residual()`](https://rdrr.io/r/stats/df.residual.html)) for
  compatibility with tools like `lmtest::coeftest()`.
- Added custom methods to indicate unsupported standard likelihood and
  residual-based generics (e.g.,
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`profile()`](https://rdrr.io/r/stats/profile.html),
  [`anova()`](https://rdrr.io/r/stats/anova.html),
  [`extractAIC()`](https://rdrr.io/r/stats/extractAIC.html),
  [`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html),
  [`rstudent()`](https://rdrr.io/r/stats/influence.measures.html)).
- Consolidated shared printing behaviors for adjustment objects
  (`adjELE`, `adjMixture`, `adjMixBayes`) using a unified base
  [`print.adjustment()`](https://postlink-group.github.io/postlink/reference/print.adjustment.md)
  method via [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html),
  reducing code duplication.
- Added an error message to clarify that the `gaussian` family in
  `glmMixture` is intended to be used only with the `identity` link for
  now.

#### Documentation & Testing

- Standardized all Roxygen manual page titles to Title Case. Implemented
  a GitHub Actions workflow (`format-titles.yml`) for documentation
  styling in continuous integration.
- Removed `VignetteBuilder: knitr` from the `DESCRIPTION` file to
  resolve a CRAN NOTE, as the extended package articles are hosted only
  via `pkgdown` for now.
- Updated the `testthat` suite to align with the newly refactored S3
  class structures and console outputs.

## **postlink 0.1.0**

CRAN release: 2026-04-15

- Initial CRAN release.
- Implements a suite of statistical tools (weighting or mixture
  modeling) for secondary analysis of linked data accounting for
  mismatch errors.
- Added support for adjusting generalized linear models, Cox
  proportional hazards models, parametric survival modeling, and
  contingency tables.
