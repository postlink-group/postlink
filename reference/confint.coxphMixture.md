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
set.seed(201)

# Simulate survival data (N = 200)
n <- 200
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.5)
true_time <- rexp(n, rate = exp(0.5 * x1 - 0.5 * x2))
cens_time <- rexp(n, rate = 0.5)

# Simulate auxiliary match scores and linkage errors
# Lower match scores correspond to a higher probability of mismatch
match_score <- runif(n, 0.5, 1.0)
is_mismatch <- rbinom(n, 1, prob = 1 - match_score)

# Induce linkage errors by shuffling covariates of mismatched records
linked_x1 <- x1
linked_x2 <- x2
mis_idx <- which(is_mismatch == 1)
shuffled_idx <- sample(mis_idx)
linked_x1[mis_idx] <- x1[shuffled_idx]
linked_x2[mis_idx] <- x2[shuffled_idx]

linked_data <- data.frame(
  time = pmin(true_time, cens_time),
  status = as.numeric(true_time <= cens_time),
  x1 = linked_x1, x2 = linked_x2,
  match_score = match_score
)

# Fit the Cox PH Mixture Model (Slawski et al., 2023)
adj <- adjMixture(linked.data = linked_data, m.formula = ~ match_score)
fit <- plcoxph(Surv(time, status) ~ x1 + x2, adjustment = adj,
               control = list(max.iter = 15))

# Extract Confidence Intervals
confint(fit)
#>                  2.5 %     97.5 %
#> x1           0.6710011  1.6021040
#> x2          -2.1429623 -0.4823973
#> (Intercept)         NA         NA
#> match_score         NA         NA
confint(fit, parm = "treatment", level = 0.90)
#>           5 % 95 %
#> treatment  NA   NA
```
