# Variance-Covariance Matrix for coxphELE Objects

Extracts the variance-covariance matrix of the main parameters (outcome
model coefficients) from a fitted `coxphELE` model.

## Usage

``` r
# S3 method for class 'coxphELE'
vcov(object, ...)
```

## Arguments

- object:

  An object of class `coxphELE`.

- ...:

  Additional arguments (currently ignored).

## Value

A matrix of the estimated covariances between the parameter estimates.

## Examples

``` r
library(survival)
set.seed(101)

# 1. Simulate data based on Vo et al. (2024) Section 3.1
n <- 200
x1 <- rnorm(n, 0, 1)
x2 <- rbinom(n, 1, 0.7)
true_hazard <- exp(0.5 * x1 - 0.5 * x2)
true_time <- rexp(n, rate = true_hazard)
cens_time <- rexp(n, rate = 0.5)
obs_time <- pmin(true_time, cens_time)
obs_status <- as.numeric(true_time <= cens_time)

# 2. Induce 15% non-informative linkage error (ELE model)
m_rate <- 0.15
mis_idx <- sample(1:n, size = round(m_rate * n))
linked_x1 <- x1
linked_x2 <- x2
shuff_idx <- sample(mis_idx)
linked_x1[mis_idx] <- x1[shuff_idx]
linked_x2[mis_idx] <- x2[shuff_idx]

linked_data <- data.frame(time = obs_time, status = obs_status,
                          x1 = linked_x1, x2 = linked_x2)

# 3. Create adjustment object and fit the model
adj <- adjELE(linked.data = linked_data, m.rate = m_rate)
fit <- plcoxph(Surv(time, status) ~ x1 + x2, adjustment = adj)

# 4. Extract the sandwich variance-covariance matrix
vmat <- vcov(fit)
print(vmat)
#>              x1           x2
#> x1  0.016220782 -0.001987996
#> x2 -0.001987996  0.083018655
```
