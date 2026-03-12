# Fit Generalized Linear Models with Linkage Error Adjustment

`plglm` fits generalized linear models (GLMs) to linked data,
incorporating adjustments for linkage error as specified in the provided
`adjustment` object. It mimics the interface of
[`glm`](https://rdrr.io/r/stats/glm.html) to ensure familiarity for
users.

## Usage

``` r
plglm(
  formula,
  family = gaussian,
  adjustment,
  subset,
  na.action,
  model = TRUE,
  x = FALSE,
  y = FALSE,
  control = list(),
  ...
)
```

## Arguments

- formula:

  A symbolic description of the model to be fitted.

- family:

  A description of the error distribution and link function to be used
  in the model. This can be a character string naming a family function,
  a family function or the result of a call to a family function.

- adjustment:

  An object inheriting from class `"adjustment"` (e.g., `adjELE`,
  `adjMixture`), or a `list` containing the necessary data and parameter
  specifications.

- subset:

  An optional vector specifying a subset of observations to be used in
  the fitting process.

- na.action:

  A function which indicates what should happen when the data contain
  NAs. The default is set by the `na.action` setting of `options`, and
  is [`na.fail`](https://rdrr.io/r/stats/na.fail.html) if that is unset.

- model:

  Logical; if `TRUE` (default), the model frame is returned.

- x, y:

  Logical; if `TRUE`, the model matrix (`x`) and response vector (`y`)
  are returned. Default is `FALSE` for both.

- control:

  A list of parameters for controlling the linkage error adjustment
  process.

- ...:

  Additional arguments passed to the underlying fitting function.

## Value

An object representing the fitted model. The specific class and
structure of the returned object depend directly on the `adjustment`
method provided:

- If `adjustment` is of class `adjELE`, returns an object of class
  [`glmELE`](https://postlink-group.github.io/postlink/reference/glmELE.md).

- If `adjustment` is of class `adjMixture`, returns an object of class
  [`glmMixture`](https://postlink-group.github.io/postlink/reference/glmMixture.md).

- If `adjustment` is of class `adjMixBayes`, returns an object of class
  [`glmMixBayes`](https://postlink-group.github.io/postlink/reference/glmMixBayes.md).

## Details

This function attempts to extract the linked data from the `adjustment`
object. It supports both reference-based storage (via `data_ref`) and
direct list components (`adjustment$data`). If the data is not present
(e.g., NULL), the function will attempt to resolve variables from the
environment of the `formula`.

It applies the standard `model.frame` processing steps (formula parsing,
subsetting, NA handling) and dispatches the resulting design matrix and
response vector to the appropriate `fitglm` method.

## See also

[`adjELE`](https://postlink-group.github.io/postlink/reference/adjELE.md),
[`adjMixture`](https://postlink-group.github.io/postlink/reference/adjMixture.md),
[`adjMixBayes`](https://postlink-group.github.io/postlink/reference/adjMixBayes.md),
[`glmELE`](https://postlink-group.github.io/postlink/reference/glmELE.md),
[`glmMixture`](https://postlink-group.github.io/postlink/reference/glmMixture.md),
[`glmMixBayes`](https://postlink-group.github.io/postlink/reference/glmMixBayes.md)

## Examples

``` r
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

# Phase 2: Estimation & Inference
# Fit a Gaussian regression model utilizing a cubic polynomial for year of birth.
fit <- plglm(
 age_at_death ~ poly(unit_yob, 3, raw = TRUE),
 family = "gaussian",
 adjustment = adj_object
)
```
