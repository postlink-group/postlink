# Model Predictions for coxphMixture Objects

Compute fitted values and predictions for the outcome model component of
the mixture. The predictions are conditional on the latent status being
a "correct match".

## Usage

``` r
# S3 method for class 'coxphMixture'
predict(
  object,
  newdata,
  type = c("lp", "risk", "expected", "survival"),
  se.fit = FALSE,
  na.action = stats::na.pass,
  reference = "strata",
  ...
)
```

## Arguments

- object:

  An object of class `coxphMixture`.

- newdata:

  Optional new data frame. If missing, predictions are for the original
  data.

- type:

  The type of prediction.

  - `"lp"`: Linear predictor (eta = X \* beta).

  - `"risk"`: Risk score (exp(eta)).

  - `"expected"`: Expected number of events (approximate).

  - `"survival"`: Survival probability at the observed times.

- se.fit:

  Logical; whether to compute standard errors (based on the
  sandwich/Louis variance).

- na.action:

  Function to handle missing values in `newdata`.

- reference:

  Reference for centering (currently ignored, defaults to uncentered).

- ...:

  Additional arguments passed to methods.

## Value

A vector or matrix of predictions, or a list containing `fit` and
`se.fit` if standard errors are requested.

## Details

When `newdata` is supplied, the function constructs the model matrix
using the terms from the original fit. Standard errors are computed
using the estimated variance-covariance matrix of the mixture model
coefficients.

For `type = "expected"` and `"survival"`, the function reconstructs the
cumulative baseline hazard step function \\\Lambda_0(t)\\ using the
Breslow estimator stored in the object and evaluates it at the time
points found in `newdata`.

## Examples

``` r
library(survival)
set.seed(205)

# Simulate auxiliary match scores and linkage errors
# Lower match scores correspond to a higher probability of mismatch
n <- 200
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.5)
true_time <- rexp(n, rate = exp(0.5 * x1 - 0.5 * x2))
cens_time <- rexp(n, rate = 0.5)
match_score <- runif(n, 0.5, 1.0)

linked_data <- data.frame(
  time = pmin(true_time, cens_time),
  status = as.numeric(true_time <= cens_time),
  x1 = x1, x2 = x2, match_score = match_score
)

mis_idx <- which(rbinom(n, 1, prob = 1 - match_score) == 1)
linked_data$x1[mis_idx] <- linked_data$x1[sample(mis_idx)]
linked_data$x2[mis_idx] <- linked_data$x2[sample(mis_idx)]

# Fit the Cox PH Mixture Model
# Note: We set `y = TRUE` to store the response for baseline hazard reconstruction
adj <- adjMixture(linked.data = linked_data, m.formula = ~ match_score)
fit <- plcoxph(Surv(time, status) ~ x1 + x2, adjustment = adj,
               y = TRUE, control = list(max.iter = 15))

# 1. Extract linear predictors for the training data
lp_train <- predict(fit, type = "lp")
head(lp_train)
#> [1]  2.2995675 -1.3674625 -0.9769177 -3.3279029  1.2543794 -2.4888121

# 2. Predict hazard ratios (risk) and expected events for a new cohort
new_cohort <- data.frame(
  time = c(1.0, 2.5, 0.5), # Required for expected/survival types
  x1 = c(0, 1.5, -1),
  x2 = c(0, 1, 1)
)

# Predict risk (exp(lp))
risk_scores <- predict(fit, newdata = new_cohort, type = "risk")
print(risk_scores)
#> [1] 1.00000000 1.27412951 0.07931427

# Predict expected number of events based on the reconstructed baseline hazard
exp_events <- predict(fit, newdata = new_cohort, type = "expected")
print(exp_events)
#> [1] 1.02610942 4.41597111 0.03660912
```
