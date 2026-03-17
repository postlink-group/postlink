# Summary method for survMixBayes models

Computes posterior summaries for the regression coefficients, mixing
weight, and component-specific distribution parameters in a fitted
`survMixBayes` model. Throughout, component 1 is interpreted as the
correct-match component and component 2 as the incorrect-match
component.

## Usage

``` r
# S3 method for class 'survMixBayes'
summary(object, probs = c(0.025, 0.5, 0.975), ...)
```

## Arguments

- object:

  An object of class `survMixBayes`.

- probs:

  Numeric vector of probabilities used to compute posterior quantiles
  for the model parameters. The default, `c(0.025, 0.5, 0.975)`, gives a
  posterior median and a 95\\ credible interval.

- ...:

  Further arguments (unused).

## Value

An object of class `summary.survMixBayes` containing posterior quantile
summaries for the regression coefficients in both mixture components,
the mixing weight, and any family-specific distribution parameters
included in the fitted model. Component 1 corresponds to the
correct-match component and component 2 to the incorrect-match
component.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(301)
n <- 150
trt <- rbinom(n, 1, 0.5)

# Simulate Weibull AFT data
true_time <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))
cens_time <- rexp(n, rate = 0.1)
true_obs_time <- pmin(true_time, cens_time)
true_status <- as.integer(true_time <= cens_time)

# Induce linkage mismatch errors in approximately 20% of records
is_mismatch <- rbinom(n, 1, 0.2)
obs_time <- true_obs_time
obs_status <- true_status
mismatch_idx <- which(is_mismatch == 1)

if (length(mismatch_idx) > 1) {
  shuffled <- sample(mismatch_idx)
  obs_time[mismatch_idx] <- obs_time[shuffled]
  obs_status[mismatch_idx] <- obs_status[shuffled]
}

linked_df <- data.frame(time = obs_time, status = obs_status, trt = trt)
adj <- adjMixBayes(linked.data = linked_df)

fit <- plsurvreg(
  survival::Surv(time, status) ~ trt,
  dist = "weibull",
  adjustment = adj,
  control = list(iterations = 200, burnin.iterations = 100, seed = 123)
)

fit_summary <- summary(fit, probs = c(0.025, 0.5, 0.975))
print(fit_summary)
} # }
```
