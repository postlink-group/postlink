# Print a glmMixture Object

Print a glmMixture Object

## Usage

``` r
# S3 method for class 'glmMixture'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `glmMixture`.

- digits:

  The number of significant digits to use.

- ...:

  Additional arguments.

## Examples

``` r
# Load the LIFE-M demo dataset
data(lifem)

# Phase 1: Adjustment Specification
# We model the correct match indicator via logistic regression using
# name commonness scores (commf, comml) and a 5% expected mismatch rate.
adj_object <- adjMixture(
 linked.data = lifem,
 m.formula = ~ commf + comml,
 m.rate = 0.05,
 safe.matches = hndlnk
)

# Phase 2: Estimation & Inference
# Fit a Gaussian regression model utilizing a cubic polynomial for year of birth.
fit <- plglm(
 age_at_death ~ poly(unit_yob, 3, raw = TRUE),
 family = "gaussian",
 adjustment = adj_object
)

print(fit)
#> 
#> Call:  plglm(formula = age_at_death ~ poly(unit_yob, 3, raw = TRUE), 
#>     family = "gaussian", adjustment = adj_object)
#> 
#> Coefficients (Outcome Model):
#>                    (Intercept)  poly(unit_yob, 3, raw = TRUE)1  
#>                          57.75                          -43.76  
#> poly(unit_yob, 3, raw = TRUE)2  poly(unit_yob, 3, raw = TRUE)3  
#>                         114.90                          -57.14  
#> 
#> Coefficients (Mismatch Model):
#> (Intercept)        commf        comml  
#>      -7.562        6.731        8.974  
#> 
#> Dispersion parameter estimate:  373.1 
#> 
```
