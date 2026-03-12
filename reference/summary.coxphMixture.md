# Summarizing Cox PH Mixture Fits

`summary` method for class `coxphMixture`. Provides a detailed summary
of the fitted model, including coefficients, hazard ratios, standard
errors, z-statistics, and p-values for both the outcome model and the
mismatch model.

## Usage

``` r
# S3 method for class 'coxphMixture'
summary(object, conf.int = 0.95, scale = 1, ...)
```

## Arguments

- object:

  An object of class `coxphMixture`.

- conf.int:

  The confidence level for the confidence intervals of the hazard
  ratios.

- scale:

  Scale factor for the standard errors (default is 1).

- ...:

  Additional arguments.

## Value

An object of class `summary.coxphMixture` containing:

- call:

  The function call.

- n:

  Total number of observations.

- nevent:

  Number of events.

- coefficients:

  Matrix of coefficients for the outcome model.

- m.coefficients:

  Matrix of coefficients for the mismatch model.

- conf.int:

  Matrix of confidence intervals for the hazard ratios.

- logtest:

  Log-likelihood information (Outcome Model).

- avgcmr:

  The average posterior probability of a correct match.

## Examples

``` r
library(survival)
set.seed(203)

# Simulate linked data with mismatch errors
# Lower match scores correspond to a higher probability of mismatch
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
linked_data$x1[mis_idx] <- linked_data$x1[sample(mis_idx)]
linked_data$x2[mis_idx] <- linked_data$x2[sample(mis_idx)]

# Fit the Cox PH Mixture Model
adj <- adjMixture(linked.data = linked_data, m.formula = ~ match_score)
fit <- plcoxph(Surv(time, status) ~ x1 + x2, adjustment = adj,
               control = list(max.iter = 15))

# Print detailed statistical summary
# Includes coefficients, hazard ratios, and the estimated average correct match rate
sum_fit <- summary(fit)
print(sum_fit)
#> 
#> Call:
#> plcoxph(formula = Surv(time, status) ~ x1 + x2, adjustment = adj, 
#>     control = list(max.iter = 15))
#> 
#> --- Outcome Model (Cox PH) ---
#>       coef exp(coef) se(coef)      z Pr(>|z|)    
#> x1  1.7950    6.0192   0.2992  6.000 1.98e-09 ***
#> x2 -1.8490    0.1574   0.4910 -3.766 0.000166 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> --- Hazard Ratios & Confidence Intervals ---
#>    exp(coef) exp(-coef) lower .95 upper .95
#> x1    6.0192     0.1661   3.34872    10.819
#> x2    0.1574     6.3537   0.06012     0.412
#> 
#> --- Mismatch Indicator Model ---
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)   -1.422      1.369  -1.039    0.299
#> match_score    1.571      1.861   0.844    0.398
#> 
#> Average Estimated Correct Match Rate: 0.4438 
#> Events: 129  / Total: 200 
#> Iterations: 15 
#> 
```
