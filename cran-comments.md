

## Submission Notes
DRAFT for the next release (to be completed by the maintainer: set the version
number in DESCRIPTION and NEWS.md and refresh the check results below).

This release replaces the Stan (rstan) backend of the Bayesian mixture models
(`glmMixBayes()`, `survregMixBayes()`, used through `adjMixBayes()`) by a Gibbs
sampler written in C++ with Rcpp / RcppArmadillo, and adds linkage information
(safe matches, a match-probability model and a prior mismatch rate) to
`adjMixBayes()`. See NEWS.md for details.

## Changes affecting the CRAN build
* The package no longer depends on `rstan`, `StanHeaders`, `rstantools`,
  `RcppParallel`, `RcppEigen`, `BH` or `label.switching`; it now links to `Rcpp`
  and `RcppArmadillo`, and `coda` is suggested.
* The Stan model files, the generated `stanExports_*` sources, `configure`,
  `configure.win` and the Windows-specific compiler flags are gone;
  `src/Makevars` and `src/Makevars.win` only link BLAS/LAPACK.
* `SystemRequirements: GNU make` was removed.
* The installed size dropped from about 6.3 Mb to about 1.9 Mb, so the three
  NOTEs of the previous submission (installed size, GNU make, `-Wa,-mbig-obj`)
  no longer occur.
* `Depends: R (>= 3.6.0)` (delayed S3 registration of the `coda::as.mcmc()`
  methods).

## Test environments
* Local: Windows 11, R 4.4.1 (Rtools 4.4)
* To be completed: GitHub Actions (macOS, Windows, Ubuntu devel / release /
  oldrel-1) and win-builder (devel)

## R CMD check results
Local `R CMD check --as-cran`: 0 errors | 0 warnings | 2 notes

* Version contains large components (0.1.2.9000): development version number,
  disappears with the release version number.
* Unable to verify current time: specific to the checking machine.

## Reverse dependencies
There are currently no reverse dependencies for this package.
