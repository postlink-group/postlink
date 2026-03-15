# Predict for survMixBayes

Returns posterior mean linear predictors for each component. Component 1
is interpreted as the correct-match component and component 2 as the
incorrect-match component.

## Usage

``` r
# S3 method for class 'survMixBayes'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  An object of class `survMixBayes`.

- newdata:

  Optional design matrix for prediction. If `NULL`, uses `object$X` if
  present.

- ...:

  Additional arguments (unused).

## Value

A list with posterior mean linear predictors for each component:
`component1` for the correct-match component and `component2` for the
incorrect-match component.

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

# Create a new design matrix for prediction
X_new <- stats::model.matrix(~ trt, data = data.frame(trt = c(0, 1)))

# Predict posterior mean linear predictors for each latent component
preds <- predict(fit, newdata = X_new)
print(preds$component1)
print(preds$component2)
} # }
```
