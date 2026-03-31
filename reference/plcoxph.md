# Fit Cox Proportional Hazards Models with Linkage Error Adjustment

`plcoxph` fits Cox proportional hazards models to linked data, adjusting
for potential mismatch errors. It serves as a wrapper around the
internal `fitcoxph` function, for compatibility with the
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) syntax.

## Usage

``` r
plcoxph(
  formula,
  adjustment,
  subset,
  na.action,
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

- model:

  Logical; if `TRUE`, the model frame is returned.

- x, y:

  Logical; if `TRUE`, the model matrix (`x`) and response (`y`) are
  returned. Defaults are `FALSE` and `FALSE`.

- control:

  A list of parameters for controlling the linkage error adjustment
  process.

- ...:

  Additional arguments passed to the internal function.

## Value

An object representing the fitted model. The specific class and
structure of the returned object depend directly on the `adjustment`
method provided:

- If `adjustment` is of class `adjELE`, returns an object of class
  [`coxphELE`](https://postlink-group.github.io/postlink/reference/coxphELE.md).

- If `adjustment` is of class `adjMixture`, returns an object of class
  [`coxphMixture`](https://postlink-group.github.io/postlink/reference/coxphMixture.md).

## See also

[`adjELE`](https://postlink-group.github.io/postlink/reference/adjELE.md),
[`adjMixture`](https://postlink-group.github.io/postlink/reference/adjMixture.md),
[`coxphELE`](https://postlink-group.github.io/postlink/reference/coxphELE.md),
[`coxphMixture`](https://postlink-group.github.io/postlink/reference/coxphMixture.md)

## Examples

``` r
library(survival)
set.seed(101)
n <- 250

# Simulate true survival data
x <- rnorm(n)
true_hazard <- exp(0.5 * x)
true_time <- rexp(n, true_hazard)
true_status <- rbinom(n, 1, 0.8)

# Induce linkage mismatch errors
match_score <- rbeta(n, 5, 1)
is_mismatch <- rbinom(n, 1, 1 - match_score)

obs_time <- true_time
obs_status <- true_status
mismatch_idx <- which(is_mismatch == 1)

# Shuffle time and status together for mismatched records
shuffled_idx <- sample(mismatch_idx)
obs_time[mismatch_idx] <- obs_time[shuffled_idx]
obs_status[mismatch_idx] <- obs_status[shuffled_idx]

linked_data <- data.frame(time = obs_time, status = obs_status, x = x, match_score)

# Specify the Adjustment Method
adj <- adjMixture(
  linked.data = linked_data,
  m.formula = ~ match_score
)

# Fit the Adjusted Cox Proportional Hazards Model
fit <- plcoxph(
  Surv(time, status) ~ x,
  adjustment = adj,
  control = list(max.iter = 50)
)
```
