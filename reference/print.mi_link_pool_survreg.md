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

dat <- data.frame(
 time = y[, "time"],
 event = y[, "event"],
 X
 )

pooled_obj <- mi_with(
  object = fit,
  data = linked_df,
  formula = survival::Surv(time, status) ~ trt
)

print(pooled_obj, digits = 4)
} # }
```
