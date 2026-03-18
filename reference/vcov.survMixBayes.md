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

shuffled <- sample(mismatch_idx)
obs_time[mismatch_idx] <- obs_time[shuffled]
obs_status[mismatch_idx] <- obs_status[shuffled]

linked_df <- data.frame(time = obs_time, status = obs_status, trt = trt)

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
