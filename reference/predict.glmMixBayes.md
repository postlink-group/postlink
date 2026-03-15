# Predictions from a glmMixBayes model

Predictions from a glmMixBayes model

## Usage

``` r
# S3 method for class 'glmMixBayes'
predict(
  object,
  newx,
  type = c("link", "response"),
  se.fit = FALSE,
  interval = c("none", "credible"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `glmMixBayes` model object.

- newx:

  A numeric matrix of new observations (n_new x K) with columns aligned
  to the design matrix `X` used for fitting.

- type:

  Either `"link"` or `"response"`, indicating the scale of predictions.

- se.fit:

  Logical; if `TRUE`, also return posterior SD of predictions.

- interval:

  Either `"none"` or `"credible"`, indicating whether to compute a
  credible interval.

- level:

  Probability level for the credible interval (default 0.95).

- ...:

  Not used.

## Value

If `se.fit = FALSE` and `interval = "none"`, a numeric vector of
predicted values. Otherwise, a matrix with columns for the fit,
(optional) `se.fit`, and (optional) credible interval bounds.

## Examples

``` r
if (FALSE) { # \dontrun{
data(lifem)

# lifem data preprocessing
# For computational efficiency in the example, we work with a subset of the lifem data.
lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
lifem_small <- rbind(
  head(subset(lifem, hndlnk == 1), 100),
  head(subset(lifem, hndlnk == 0), 20)
)

x <- cbind(1, poly(lifem_small$unit_yob, 3, raw = TRUE))
y <- lifem_small$age_at_death

adj <- adjMixBayes(
  linked.data = lifem_small,
  priors = list(theta = "beta(2, 2)")
)

fit <- plglm(
  age_at_death ~ poly(unit_yob, 3, raw = TRUE),
  family = "gaussian",
  adjustment = adj,
  control = list(
    iterations = 200,
    burnin.iterations = 100,
    seed = 123
  )
)

newx <- cbind(1, poly(c(0.2, 0.5, 0.8), 3, raw = TRUE))
predict(fit, newx = newx, type = "response")
} # }
```
