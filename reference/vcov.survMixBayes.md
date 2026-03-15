# Posterior covariance matrix for survMixBayes coefficients

Returns the empirical posterior covariance matrix of the regression
coefficients for component 1 of a fitted `survMixBayes` model. In this
package, component 1 is interpreted as the correct-match component.

## Usage

``` r
# S3 method for class 'survMixBayes'
vcov(object, ...)
```

## Arguments

- object:

  A `survMixBayes` model object.

- ...:

  Further arguments (unused).

## Value

Posterior covariance matrix of the regression coefficients for component
1, interpreted as the correct-match component.

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

# Extract the empirical posterior covariance matrix for component 1
vcov_mat <- vcov(fit)
print(vcov_mat)
} # }
```
