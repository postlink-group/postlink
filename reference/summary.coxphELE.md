# Summary of a coxphELE Object

Produces a summary of the fitted CoxPH model with linkage error
adjustment, including coefficient estimates, hazard ratios, standard
errors, z-statistics, and p-values.

## Usage

``` r
# S3 method for class 'coxphELE'
summary(object, conf.int = 0.95, ...)
```

## Arguments

- object:

  An object of class `coxphELE`.

- conf.int:

  The confidence level for the confidence intervals (default is 0.95).

- ...:

  Additional arguments (currently ignored).

## Value

An object of class `summary.coxphELE` containing:

- call:

  The matched call.

- coefficients:

  A matrix with columns for coefficients, hazard ratios (exp(coef)),
  standard errors, z-values, and p-values.

- conf.int:

  A matrix of confidence intervals for the hazard ratios.

- n:

  The number of observations.

- nevent:

  The number of events.

## Examples

``` r
library(survival)
set.seed(103)

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

# Generate and print the detailed statistical summary
sum_fit <- summary(fit)
print(sum_fit)
#> Call:
#> plcoxph(formula = Surv(time, status) ~ x1 + x2, adjustment = adj)
#> 
#> n=200, number of events=86
#> 
#>        coef exp(coef) se(coef)      z Pr(>|z|)
#> x1 -0.05155   0.94976  0.13096 -0.394    0.694
#> x2  0.35668   1.42859  0.39408  0.905    0.365
#> 
#> Confidence Intervals for Hazard Ratios:
#>    exp(coef) lower 0.95 upper 0.95
#> x1    0.9498     0.7348      1.228
#> x2    1.4286     0.6599      3.093
#> 
```
