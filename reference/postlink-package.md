# postlink: Post-Linkage Data Analysis Under Record Linkage Errors

The `postlink` package provides a unified suite of statistical tools
designed to rigorously account for record linkage errors in downstream
modeling.

Record linkage is rarely perfect. When datasets are merged using noisy
or non-unique identifiers, mismatches (false links) are inadvertently
introduced. Ignoring these errors acts as a contaminant in regression
analysis, typically leading to significantly attenuated estimates and
biased statistical inference. `postlink` equips researchers with
methodologies to propagate linkage uncertainty into their models,
specifically accommodating "secondary analysis" workflows where direct
access to the primary, unlinked files is restricted.

## Details

**A Two-Phase Workflow**

The package is built on a modular, object-oriented S3 architecture that
decouples the specification of linkage error from the substantive
statistical modeling. This provides a familiar, standard formula-based
modeling interface.

**Phase 1: Adjustment Specification**

The analyst encapsulates the linked data and the chosen error-adjustment
methodology by using one of the `adj*` constructor functions. These
functions validate the data and return a structured adjustment object:

- [`adjELE()`](https://postlink-group.github.io/postlink/reference/adjELE.md):
  Specifies the Exchangeable Linkage Error (ELE) model, which corrects
  for bias using known or audited mismatch rates.

- [`adjMixture()`](https://postlink-group.github.io/postlink/reference/adjMixture.md):
  Specifies a frequentist mixture model approach that treats match
  status as a latent variable, estimating error rates directly from data
  (e.g., using match scores) via the EM algorithm.

- [`adjMixBayes()`](https://postlink-group.github.io/postlink/reference/adjMixBayes.md):
  Specifies a Bayesian mixture model approach, enabling parameter
  estimation and multiple imputation of latent match statuses using
  Stan.

**Phase 2: Estimation & Inference**

The adjustment object is subsequently passed to a standard modeling
wrapper, integrating the error correction into the familiar R modeling
syntax:

- [`plglm()`](https://postlink-group.github.io/postlink/reference/plglm.md):
  Generalized Linear Models (linear, logistic, Poisson, Gamma).

- [`plcoxph()`](https://postlink-group.github.io/postlink/reference/plcoxph.md):
  Cox Proportional Hazards models.

- [`plctable()`](https://postlink-group.github.io/postlink/reference/plctable.md):
  Contingency table analyses.

- [`plsurvreg()`](https://postlink-group.github.io/postlink/reference/plsurvreg.md):
  Parametric survival models.

As a note, estimation and inference supported for each type of
adjustment object vary. Refer to the `adj*` documentation for models
currently supported.

Post-estimation tools function as expected in standard R workflows
(e.g., [`summary()`](https://rdrr.io/r/base/summary.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html), and
[`confint()`](https://rdrr.io/r/stats/confint.html)). These methods
specially are derived to account for the additional steps introduced by
the linkage error adjustment.

**Note:** While the two-phase workflow is recommended for standard
analyses, the package's architecture intentionally isolates the core
logic of each method into independent internal routines. Advanced users,
developers, or researchers running large-scale simulations can bypass
the wrapper functions and formula-parsing overhead entirely. By
supplying pre-computed design matrices and response vectors, the
underlying computational functions can be directly used if preferred
(e.g.,
[`coxphELE()`](https://postlink-group.github.io/postlink/reference/coxphELE.md),
[`glmMixBayes()`](https://postlink-group.github.io/postlink/reference/glmMixBayes.md),
[`survregMixBayes()`](https://postlink-group.github.io/postlink/reference/survregMixBayes.md)).

## References

Chambers, R. (2009). Regression analysis of probability-linked data.
*Official Statistics Research Series*, 4, 1-15.

Gutman, R., Sammartino, C., Green, T., & Montague, B. (2016). Error
adjustments for file linking methods using encrypted unique client
identifier (eUCI) with application to recently released prisoners who
are HIV+. *Statistics in Medicine*, 35(1), 115–129.
[doi:10.1002/sim.6586](https://doi.org/10.1002/sim.6586)

Slawski, M., West, B., Bukke, P., Wang, Z., Diao, G., & Ben-David, E.
(2023). A General Framework for Regression with Mismatched Data Based on
Mixture Modeling. *arXiv preprint arXiv:2306.00909*.

Vo, T. H., Garès, V., Zhang, L. C., Happe, A., Oger, E., Paquelet, S., &
Chauvet, G. (2024). Cox regression with linked data. *Statistics in
Medicine*, 43(2), 296-314.
[doi:10.1002/sim.9960](https://doi.org/10.1002/sim.9960)

## See also

Useful links:

- <https://postlink-group.github.io/postlink/>

- Report bugs at <https://github.com/postlink-group/postlink/issues>

## Author

**Maintainer**: Priyanjali Bukke <postlink.group@gmail.com>

Authors:

- Gauri Kamat

- Jiahao Cui

- Roee Gutman

- Martin Slawski

Other contributors:

- Zhenbang Wang \[contributor\]

- Brady T. West \[contributor\]

- Emanuel Ben-David \[contributor\]

- Guoqing Diao \[contributor\]
