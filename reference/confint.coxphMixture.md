# Confidence Intervals for coxphMixture Objects

Computes Wald confidence intervals for the coefficients of the outcome
model and the mismatch indicator model.

## Usage

``` r
# S3 method for class 'coxphMixture'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `coxphMixture`.

- parm:

  A specification of which parameters are to be given confidence
  intervals, either a vector of numbers or names. If missing, all
  parameters are considered.

- level:

  The confidence level required (default is 0.95).

- ...:

  Additional arguments (currently ignored).

## Value

A matrix (or vector) with lower and upper confidence limits for each
parameter.

## Details

The intervals are calculated based on the variance-covariance matrix
returned by
[`vcov.coxphMixture`](https://postlink-group.github.io/postlink/reference/vcov.coxphMixture.md),
using the standard normal approximation: `Estimate +/- z_crit * SE`.

## Examples

``` r
library(survival)
set.seed(202)

# Simulate linked data with heterogeneous mismatch errors
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
if (length(mis_idx) > 1) {
  linked_data$x1[mis_idx] <- linked_data$x1[sample(mis_idx)]
  linked_data$x2[mis_idx] <- linked_data$x2[sample(mis_idx)]
}

# Fit the Cox PH Mixture Model
adj <- adjMixture(linked.data = linked_data, m.formula = ~ match_score)
fit <- plcoxph(Surv(time, status) ~ x1 + x2, adjustment = adj,
               control = list(max.iter = 15))

# Extract 95% Confidence Intervals for all parameters
confint(fit)
#>                  2.5 %    97.5 %
#> x1           0.8636149 2.1058154
#> x2          -1.6027959 0.4152174
#> (Intercept)         NA        NA
#> match_score         NA        NA

# Extract 90% Confidence Intervals for a specific outcome parameter
confint(fit, parm = "x1", level = 0.90)
#>          5 %     95 %
#> x1 0.9634713 2.005959
```
