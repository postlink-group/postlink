# **postlink**

The `postlink` package is dedicated to providing a unified suite of
statistical tools designed to rigorously account for record linkage
errors in post-linkage data analysis.

Record linkage is often error-prone, particularly when identifiers used
for matching records are noisy or non-unique. Mismatches (false matches)
act as a contaminant in the linked data, typically leading to attenuated
estimates in downstream analysis. The `postlink` R package currently
supports three statistical frameworks to account for potential mismatch
errors during downstream regression modeling:

- **Weighting**: Solves adjusted estimating equations using
  block-specific proportions of correct links. It assumes an
  Exchangeable Linkage Error (ELE) model (Chambers, 2009; Chambers et
  al., 2023; Vo et al., 2024).
- **Mixture Modeling**: Uses an Expectation-Maximization (EM) algorithm
  and infers the latent correct match status using a two-component
  mixture model (Slawski et al., 2025).
- **Bayesian Mixture Modeling**: Performs posterior inference via data
  augmentation, alternating between imputing the latent match status and
  updating model parameters to propagate linkage uncertainty (Gutman et
  al., 2016).

The `postlink` package currently focuses on methods for secondary
analysis, where the individual files that were linked are not
accessible. For the primary analysis setting, when individual files are
accessible, methods that perform record linkage and analysis jointly
with direct propagation of uncertainty would be more suitable.

The long-term goal of `postlink` is to extend support for a wide array
of linkage and post-linkage analysis scenarios.

## **Package Design**

The package is built on a modular, object-oriented S3 architecture that
decouples the specification of linkage error from the substantive
statistical modeling. This provides a familiar, standard formula based
modeling interface.

**Phase 1**: Adjustment Specification

First, we define the linked data and the chosen adjustment methodology
using a constructor function. These constructors validate the data and
return a lightweight S3 adjustment object.

- [`adjELE()`](https://postlink-group.github.io/postlink/reference/adjELE.md):
  Specifies the Exchangeable Linkage Error (ELE) model, using known or
  audited mismatch rates.
- [`adjMixture()`](https://postlink-group.github.io/postlink/reference/adjMixture.md):
  Specifies a frequentist mixture model approach that treats match
  status as a latent variable, estimating linkage error rates directly
  from data (e.g., using an overall mismatch rate, safe matches,
  predictors of match status, or predicted correct match probabilities).
  If no record linkage information is available, a constant mismatch
  rate is assumed.
- [`adjMixBayes()`](https://postlink-group.github.io/postlink/reference/adjMixBayes.md):
  Specifies a Bayesian mixture model approach, with default or
  informative priors, enabling parameter estimation and multiple
  imputation of latent match statuses using Stan.

**Phase 2**: Estimation & Inference

The adjustment object is subsequently passed to a standard modeling
wrapper, integrating the linkage error correction into the familiar R
modeling syntax:

- [`plglm()`](https://postlink-group.github.io/postlink/reference/plglm.md)
  for generalized linear models (linear, logistic, Poisson, Gamma).
- [`plcoxph()`](https://postlink-group.github.io/postlink/reference/plcoxph.md)
  for Cox proportional hazards regression.
- [`plsurvreg()`](https://postlink-group.github.io/postlink/reference/plsurvreg.md)
  for parametric survival models.
- [`plctable()`](https://postlink-group.github.io/postlink/reference/plctable.md)
  for contingency table analysis.

Estimation and inference supported for each type of adjustment object
vary. Please refer to the `adj*` or `pl*` documentation for models
currently supported.

Standard R workflows (e.g., summary(), predict(), vcov(), and confint())
can be used for display and processing of results. These methods
specially are derived to consider the additional steps introduced by the
linkage error adjustment.

**Note**: While the two-phase workflow is recommended for standard
analyses, the package’s architecture isolates the core logic of each
method into independent internal routines. If preferred, the underlying
computational functions can be used directly by supplying pre-computed
design matrices and response vectors (e.g.,
[`coxphELE()`](https://postlink-group.github.io/postlink/reference/coxphELE.md),
[`glmMixture()`](https://postlink-group.github.io/postlink/reference/glmMixture.md),
[`survregMixBayes()`](https://postlink-group.github.io/postlink/reference/survregMixBayes.md))

## **Installation**

The development version of `postlink` can be installed from GitHub or
locally:

``` r
# Using devtools:
# install.packages("devtools")
devtools::install_github("postlink-group/postlink")

# Or, using pak:
# install.packages("pak")
pak::pkg_install("postlink-group/postlink")
```

## **System Requirements**

Because `postlink` includes Bayesian mixture models powered by `rstan`,
installing the development version of the package from source requires a
working C++ toolchain to compile the underlying models.

Depending on your operating system, please ensure you have the following
installed:

- **Windows:** Install the version of
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/) that matches
  your current R version.
- **macOS:** Install the Xcode Command Line Tools by opening your
  terminal and running: `xcode-select --install`
- **Linux:** Ensure you have the standard C++ compiler and R development
  tools installed (e.g., `sudo apt-get install r-base-dev` on
  Ubuntu/Debian).

*(Note: Once `postlink` is officially released on CRAN, Windows and
macOS users will be able to download pre-compiled binaries, bypassing
this requirement).*

## **Quick Start**

Below is a brief example illustrating the typical workflow using
`postlink`.

We analyze the relationship between age at death and year of birth
linked using historical records from the LIFE-M project. The linked
dataset contains a mix of hand-linked records (assumed correct) and
purely machine-linked records subject to an approximate 5% mismatch
rate. Instead of fitting the standard glm model ignoring the mismatch
errors, we use the `postlink` to adjust for potential mismatches using
the entire linked dataset.

``` r
library(postlink)

# Load the LIFE-M demo dataset
data(lifem)

# Phase 1: Adjustment Specification
# We model the correct match indicator via logistic regression using 
# name commonness scores (commf, comml) and a 5% expected mismatch rate.
adj_object <- adjMixture(
  linked.data = lifem,
  m.formula = ~ commf + comml,
  m.rate = 0.05,
  safe.matches = hndlnk
)

# Check specified adjustment
print(adj_object)

# Phase 2: Estimation & Inference
# Fit a Gaussian regression model utilizing a cubic polynomial for year of birth.
fit <- plglm(
  age_at_death ~ poly(unit_yob, 3, raw = TRUE),
  family = "gaussian",
  adjustment = adj_object
)

# View model results
summary(fit)
confint(fit)
```

## **References**

Chambers, R. (2009). Regression analysis of probability-linked data. ,
4, 1-15.

Chambers, R. L., Fabrizi, E., Ranalli, M. G., Salvati, N., & Wang, S.
(2023). Robust regression using probabilistically linked data. , 15(2),
e1596.

Gutman, R., Sammartino, C., Green, T., & Montague, B. (2016). Error
adjustments for file linking methods using encrypted unique client
identifier (eUCI) with application to recently released prisoners who
are HIV+. , 35(1), 115–129.

Slawski, M., West, B. T., Bukke, P., Wang, Z., Diao, G., & Ben-David, E.
(2025). A general framework for regression with mismatched data based on
mixture modelling. , 188(3), 896-919.

Vo, T. H., Garès, V., Zhang, L. C., Happe, A., Oger, E., Paquelet, S., &
Chauvet, G. (2024). Cox regression with linked data. , 43(2), 296-314.
