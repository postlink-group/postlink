# Print pooled Cox regression results

Print pooled Cox regression results

## Usage

``` r
# S3 method for class 'mi_link_pool_survreg'
print(x, digits = max(3L, getOption("digits") - 2L), ...)
```

## Arguments

- x:

  An object of class `mi_link_pool_survreg`, typically returned by
  [`mi_with()`](https://postlink-group.github.io/postlink/reference/mi_with.md)
  for a `survMixBayes` fit.

- digits:

  the number of significant digits to print.

- ...:

  further arguments (unused).

## Value

The input `x`, invisibly.

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

pooled_obj <- mi_with(
  object = fit,
  data = linked_df,
  formula = survival::Surv(time, status) ~ trt
)

print(pooled_obj, digits = 4)
} # }
```
