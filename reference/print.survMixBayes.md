# Print a survMixBayes model object

Prints the model call and posterior mean regression coefficients for the
first mixture component of the fitted survival model. In this package,
component 1 is interpreted as the correct-match component and component
2 as the incorrect-match component.

## Usage

``` r
# S3 method for class 'survMixBayes'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `survMixBayes`.

- digits:

  Minimum number of significant digits to show.

- ...:

  Further arguments (unused).

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

print(fit)
} # }
```
