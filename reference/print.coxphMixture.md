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

## Value

Invisibly returns the input object `x`.

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

# Explicitly call the print method
print(fit)
#> Call:
#> plcoxph(formula = Surv(time, status) ~ x1 + x2, adjustment = adj, 
#>     control = list(max.iter = 15))
#> 
#> Outcome Model Coefficients:
#>     x1     x2 
#>  1.137 -1.313 
#> 
#> Mismatch Model Coefficients:
#> (Intercept) match_score 
#>      -5.687       7.128 
#> 
#> Likelihood ratio test (model=outcome) not available due to pseudo-likelihood.
#> n= 200 , number of events= 113 
```
