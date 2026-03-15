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
n <- 200

# 1. Simulate covariates
age_centered <- rnorm(n, 0, 5)
treatment <- rbinom(n, 1, 0.5)

true_time <- rexp(n, rate = exp(0.05 * age_centered - 0.5 * treatment))
cens_time <- rexp(n, rate = 0.2)
time <- pmin(true_time, cens_time)
status <- as.numeric(true_time <= cens_time)

# 2. Simulate mismatch errors based on a matching score
# Lower match scores correspond to a higher probability of mismatch
match_score <- runif(n, 0.5, 1.0)
is_mismatch <- rbinom(n, 1, prob = 1 - match_score)
mis_idx <- which(is_mismatch == 1)

linked_age <- age_centered
linked_trt <- treatment

if(length(mis_idx) > 1){
 shuffled <- sample(mis_idx)
 linked_age[mis_idx] <- age_centered[shuffled]
 linked_trt[mis_idx] <- treatment[shuffled]
}

linked_data <- data.frame(time = time, status = status,
                          age = linked_age, treatment = linked_trt,
                          match_score = match_score)

# 3. Fit the Cox PH Mixture Model
adj <- adjMixture(linked.data = linked_data, m.formula = ~ match_score)
fit <- plcoxph(Surv(time, status) ~ age + treatment, adjustment = adj,
               control = list(max.iter = 15))

# 4. Extract Confidence Intervals
confint(fit)
#>                   2.5 %      97.5 %
#> age         -0.02919245  0.11861452
#> treatment   -1.57865485 -0.07102464
#> (Intercept)          NA          NA
#> match_score          NA          NA
confint(fit, parm = "treatment", level = 0.90)
#>                 5 %       95 %
#> treatment -1.457461 -0.1922182
```
