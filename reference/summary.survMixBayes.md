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

# Component 1 represents correct links (signal),
# and component 2 represents incorrect links (noise).
Z_true <- 2 - rbinom(n, 1, 0.8)

time1 <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))  # Correct links
time2 <- rweibull(n, shape = 1.2, scale = exp(1 + 0.2 * trt))  # Incorrect links
obs_time <- ifelse(Z_true == 1, time1, time2)

cens_time <- rexp(n, rate = 0.1)
status <- as.integer(obs_time <= cens_time)
obs_time <- pmin(obs_time, cens_time)

linked_df <- data.frame(time = obs_time, status = status, trt = trt)

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
