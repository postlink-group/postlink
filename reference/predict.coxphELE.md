# Predictions from a coxphELE Object

Computes linear predictors or risk scores from a fitted `coxphELE`
model. If `newdata` is provided, the method attempts to construct the
model matrix from the call stored in the object.

## Usage

``` r
# S3 method for class 'coxphELE'
predict(object, newdata = NULL, type = c("lp", "risk"), ...)
```

## Arguments

- object:

  An object of class `coxphELE`.

- newdata:

  Optional new data frame to obtain predictions for. If omitted, the
  linear predictors of the original data are returned.

- type:

  The type of prediction required. Type `"lp"` (default) returns the
  linear predictor. Type `"risk"` returns the hazard ratio, `exp(lp)`.

- ...:

  Additional arguments (currently ignored).

## Value

A numeric vector of predictions.

## Details

If `newdata` is supplied, the function attempts to retrieve the formula
from `object$call` or `object$formula`. If the original model was fitted
using the internal engine directly without a wrapper that stores the
call/formula, prediction on new data will fail.

## Examples

``` r
library(survival)
set.seed(105)

# Simulate linked data subject to mismatch error
n <- 200
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.7)
true_time <- rexp(n, rate = exp(0.5 * x1 - 0.5 * x2))
cens_time <- rexp(n, rate = 0.5)

# Induce 15% linkage error
mis_idx <- sample(1:n, size = 0.15 * n)
x1[mis_idx] <- x1[sample(mis_idx)]
x2[mis_idx] <- x2[sample(mis_idx)]

# Linked data
linked_data <- data.frame(
  time = pmin(true_time, cens_time),
  status = as.numeric(true_time <= cens_time),
  x1 = x1, x2 = x2
)

# Fit the adjusted Cox PH model
adj <- adjELE(linked.data = linked_data, m.rate = 0.15)
fit <- plcoxph(Surv(time, status) ~ x1 + x2, adjustment = adj)

# 1. Extract linear predictors for the original data
lp_train <- predict(fit, type = "lp")
head(lp_train)
#> [1] -0.007068736 -0.042136465 -0.026928050 -0.043077001 -0.026705796
#> [6] -0.015142319

# 2. Predict hazard ratios (risk) for a new cohort
new_cohort <- data.frame(x1 = c(0, 1.5, -1), x2 = c(0, 1, 1))
risk_scores <- predict(fit, newdata = new_cohort, type = "risk")
print(risk_scores)
#> [1] 1.0000000 0.9551026 0.9889093
```
