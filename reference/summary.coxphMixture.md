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

# Print detailed statistical summary
sum_fit <- summary(fit)
print(sum_fit)
#> 
#> Call:
#> plcoxph(formula = Surv(time, status) ~ x1 + x2, adjustment = adj, 
#>     control = list(max.iter = 15))
#> 
#> --- Outcome Model (Cox PH) ---
#>       coef exp(coef) se(coef)      z Pr(>|z|)    
#> x1  1.1366    3.1160   0.2375  4.785 1.71e-06 ***
#> x2 -1.3127    0.2691   0.4236 -3.099  0.00194 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> --- Hazard Ratios & Confidence Intervals ---
#>    exp(coef) exp(-coef) lower .95 upper .95
#> x1    3.1160     0.3209    1.9562    4.9635
#> x2    0.2691     3.7161    0.1173    0.6173
#> 
#> --- Mismatch Indicator Model ---
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)   -5.687      2.337  -2.434   0.0149 *
#> match_score    7.128      2.977   2.395   0.0166 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Average Estimated Correct Match Rate: 0.4527 
#> Events: 113  / Total: 200 
#> Iterations: 15 
#> 
```
