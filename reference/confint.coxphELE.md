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
n <- 200

# 1. Simulate covariates
age_centered <- rnorm(n, 0, 5)
treatment <- rbinom(n, 1, 0.5)

# 2. Simulate true survival times
true_time <- rexp(n, rate = exp(0.05 * age_centered - 0.6 * treatment))
cens_time <- rexp(n, rate = 0.2)
time <- pmin(true_time, cens_time)
status <- as.numeric(true_time <= cens_time)

# 3. Induce 15% Exchangeable Linkage Error (ELE)
mis_idx <- sample(1:n, size = floor(0.15 * n))
linked_age <- age_centered
linked_trt <- treatment

 # False links drawn uniformly from the target population
 false_link_idx <- sample(1:n, size = length(mis_idx), replace = TRUE)
 linked_age[mis_idx] <- age_centered[false_link_idx]
 linked_trt[mis_idx] <- treatment[false_link_idx]

linked_data <- data.frame(time = time, status = status,
                          age = linked_age, treatment = linked_trt)

# 4. Fit the adjusted Cox PH model
adj <- adjELE(linked.data = linked_data, m.rate = 0.15)
fit <- plcoxph(Surv(time, status) ~ age + treatment, adjustment = adj)

# 5. Compute confidence intervals
confint(fit) # 95% CI for all coefficients
#>                2.5 %     97.5 %
#> age       -0.0853233 0.06064571
#> treatment -0.6688953 0.89980481
confint(fit, parm = "treatment", level = 0.90) # 90% CI for a specific parameter
#>                  5 %      95 %
#> treatment -0.5427925 0.7737021
```
