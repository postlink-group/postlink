# Print a coxphMixture Object

Print a coxphMixture Object

## Usage

``` r
# S3 method for class 'coxphMixture'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `coxphMixture`.

- digits:

  The number of significant digits to use.

- ...:

  Additional arguments.

## Examples

``` r
library(survival)
set.seed(204)

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
}

# Fit the Cox PH Mixture Model
adj <- adjMixture(linked.data = linked_data, m.formula = ~ match_score)
fit <- plcoxph(Surv(time, status) ~ x1 + x2, adjustment = adj,
               control = list(max.iter = 15))

# Explicitly call the print method
print(fit)
#> Call:
#> plcoxph(formula = Surv(time, status) ~ x1 + x2, adjustment = adj, 
#>     control = list(max.iter = 15))
#> 
#> Outcome Model Coefficients:
#>      x1      x2 
#>  1.6929 -0.2021 
#> 
#> Mismatch Model Coefficients:
#> (Intercept) match_score 
#>       3.543      -5.091 
#> 
#> Likelihood ratio test (model=outcome) not available due to pseudo-likelihood.
#> n= 200 , number of events= 129 
```
