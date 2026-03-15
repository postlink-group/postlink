# Pool regression fits across posterior draws of correct-match classifications

Use posterior draws of the latent match indicators from
[`survregMixBayes()`](https://postlink-group.github.io/postlink/reference/survregMixBayes.md)
to repeatedly identify which records are treated as correct matches,
refit a Cox proportional hazards model on those records, and pool the
resulting estimates using multiple-imputation pooling rules.

Each retained posterior draw defines one subset of records classified as
correct matches. The function fits the specified
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) model
to that subset, extracts the estimated coefficients and covariance
matrix, and combines the results across draws using Rubin's rules.

## Usage

``` r
# S3 method for class 'survMixBayes'
mi_with(
  object,
  data,
  formula,
  min_n = NULL,
  quietly = TRUE,
  ties = "efron",
  ...
)
```

## Arguments

- object:

  A `survMixBayes` model object containing posterior draws of the latent
  match indicators.

- data:

  A data.frame with all candidate records in the same row order as used
  in the model.

- formula:

  Model formula for refitting on each draw (required), typically of the
  form `survival::Surv(time, event) ~ ...`.

- min_n:

  Minimum number of records required to fit the model for a given
  posterior draw. The default is `p + 2`, where `p` is the number of
  non-intercept columns in the model matrix.

- quietly:

  If `TRUE`, draws that lead to fitting errors are skipped without
  printing the full error message.

- ties:

  Method for handling tied event times in
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).
  Default is `"efron"`.

- ...:

  Additional arguments passed to
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).

## Value

An object of class `c("mi_link_pool_survreg", "mi_link_pool")`
containing pooled coefficient estimates, standard errors, confidence
intervals, and related summary information.

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

print(pooled_obj)
} # }
```
