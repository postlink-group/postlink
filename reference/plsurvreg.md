# Fit Parametric Survival Models with Linkage Error Adjustment

`plsurvreg` fits parametric survival models (Accelerated Failure Time
models) to linked data. It serves as a wrapper for the internal engines,
compatible with
[`survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
specifications.

## Usage

``` r
plsurvreg(
  formula,
  adjustment,
  subset,
  na.action,
  dist = "weibull",
  model = TRUE,
  x = FALSE,
  y = FALSE,
  control = list(),
  ...
)
```

## Arguments

- formula:

  A formula object, with the response on the left of a ~ operator, and
  the terms on the right. The response must be a survival object as
  returned by the [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)
  function.

- adjustment:

  An object inheriting from class `"adjustment"`, or a `list` containing
  the necessary parameter specifications.

- subset:

  An optional vector specifying a subset of observations.

- na.action:

  A function for handling NAs.

- dist:

  Character string specifying the survival distribution (currently it
  must be "weibull" or "gamma").

- model:

  Logical; if `TRUE`, the model frame is returned.

- x, y:

  Logical; if `TRUE`, the model matrix (`x`) and response (`y`) are
  returned.

- control:

  A list of control parameters.

- ...:

  Additional arguments passed to the internal engine.

## Value

A fitted survival model object.

## See also

[`survreg`](https://rdrr.io/pkg/survival/man/survreg.html)

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
set.seed(202)
n <- 200

# Simulate Weibull AFT data
trt <- rbinom(n, 1, 0.5)
true_time <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))
cens_time <- rexp(n, 0.05)
true_obs_time <- pmin(true_time, cens_time)
true_status <- as.numeric(true_time <= cens_time)

# Induce linkage mismatch errors
is_mismatch <- rbinom(n, 1, 0.2) # ~20% overall mismatch rate
obs_time <- true_obs_time
obs_status <- true_status
mismatch_idx <- which(is_mismatch == 1)

if(length(mismatch_idx) > 1) {
  shuffled <- sample(mismatch_idx)
  obs_time[mismatch_idx] <- obs_time[shuffled]
  obs_status[mismatch_idx] <- obs_status[shuffled]
}
linked_df <- data.frame(time = obs_time, status = obs_status, trt)

# Specify the Bayesian Mixture Adjustment
adj <- adjMixBayes(linked.data = linked_df)

# Fit the Adjusted Parametric Survival Model
fit <- plsurvreg(
  Surv(time, status) ~ trt,
  dist = "weibull",
  adjustment = adj,
  control = list(iterations = 2000, burnin.iterations = 500)
)
} # }
```
