# Confidence Intervals for coxphELE Objects

Computes Wald confidence intervals for the model coefficients.

## Usage

``` r
# S3 method for class 'coxphELE'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `coxphELE`.

- parm:

  A specification of which parameters are to be given confidence
  intervals, either a vector of numbers or names. If missing, all
  parameters are considered.

- level:

  The confidence level required (default 0.95).

- ...:

  Additional arguments (currently ignored).

## Value

A matrix (or vector) with lower and upper confidence limits for each
parameter.

## Examples

``` r
library(survival)
set.seed(104)

# Simulate linked data subject to mismatch error
n <- 200
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.7)
true_time <- rexp(n, rate = exp(0.5 * x1 - 0.5 * x2))
cens_time <- rexp(n, rate = 0.5)

linked_data <- data.frame(
  time = pmin(true_time, cens_time),
  status = as.numeric(true_time <= cens_time),
  x1 = x1, x2 = x2
)

# Induce 15% linkage error
mis_idx <- sample(1:n, size = 0.15 * n)
linked_data$x1[mis_idx] <- linked_data$x1[sample(mis_idx)]
linked_data$x2[mis_idx] <- linked_data$x2[sample(mis_idx)]

# Fit the adjusted Cox PH model
adj <- adjELE(linked.data = linked_data, m.rate = 0.15)
fit <- plcoxph(Surv(time, status) ~ x1 + x2, adjustment = adj)

# Compute 95% confidence intervals for all coefficients
confint(fit)
#>         2.5 %    97.5 %
#> x1 -0.3822596 0.1305227
#> x2 -0.3074440 1.1396112

# Compute 90% confidence intervals for a specific parameter
confint(fit, parm = "x1", level = 0.90)
#>           5 %       95 %
#> x1 -0.3410387 0.08930178
```
