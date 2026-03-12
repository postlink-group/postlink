# Contributing to postlink

Thank you for visiting this page!

This page is currently under development, but in the meantime, if you
found a bug, have a feature request, or would be interested in
contributing to the package, please open an issue or email
<postlink.group@gmail.com> to discuss the proposed changes. We look
forward to your feedback and collaborating with you.

## The `R/` Directory Layout

Files are grouped into four main categories based on their role in the
package architecture:

1.  Foundation and Utilities:

- `00_utils.R`: Shared internal helper functions used across various
  methods.
- `postlink-package.R`: The core package-level documentation and
  namespace directives.
- `stanmodels.R`: The interface for pre-compiled Stan/C++ models used in
  Bayesian adjustments.

2.  Phase 1: Adjustment Constructors (The `01_` Prefix): These files
    contain the user-facing constructor functions that define the linked
    data and specify the adjustment method.

- `01_adjELE.R`, `01_adjMixture.R`, `01_adjMixBayes.R`: The core
  constructor functions.
- `01_adjELE_methods.R`, etc.: The S3 methods (like `print`, `summary`)
  specifically for these adjustment objects.

4.  Phase 2: Wrapper Models (The `02_` Prefix): These files contain the
    user-facing modeling functions that parse formulas and dispatch the
    appropriate internal routine based on the provided adjustment
    object.

- `02_plglm.R`, `02_plcoxph.R`, `02_plctable.R`, `02_plsurvreg.R`

4.  Internal Routines for the Adjustment Methods: The internal
    computational routines are completely decoupled from the wrappers.
    They are named using a `[method]_[model]` syntax.

- ``` R
  ELE: `ele_glm.R`, `ele_coxph.R`, and their associated `*_methods.R` files.
  ```

- ``` R
  Mixture: `mixture_glm.R`, `mixture_coxph.R`, `mixture_ctable.R`, `mixture_helpers.R`, etc. 
  ```

- ``` R
  Bayesian Mixture: `mixbayes_glm.R`, `mixbayes_survreg.R`, `mixbayes_helpers.R`, etc.
  ```

5.  Datasets:

- `lifem.R`, `brfss.R`: Roxygen documentation for the included demo
  datasets.
