# Print a coxphELE Object

Prints a short summary of the fitted CoxPH model with linkage error
adjustment, displaying the call (if available) and the estimated
coefficients.

## Usage

``` r
# S3 method for class 'coxphELE'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `coxphELE`.

- digits:

  The number of significant digits to use when printing.

- ...:

  Additional arguments passed to
  [`print.default`](https://rdrr.io/r/base/print.default.html).

## Value

Invisibly returns the input object `x`.

## Examples

``` r
library(survival)
set.seed(102)

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

# Print the basic model output
print(fit)
#> Call:
#> plcoxph(formula = Surv(time, status) ~ x1 + x2, adjustment = adj)
#> 
#> Coefficients:
#>     x1      x2  
#> 0.3609  0.1438  
#> 
#> n= 200, number of events= 85
```
