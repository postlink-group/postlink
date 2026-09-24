# Summarizing GLM Mixture Fits

`summary` method for class `glmMixture`.

## Usage

``` r
# S3 method for class 'glmMixture'
summary(object, dispersion = NULL, ...)
```

## Arguments

- object:

  An object of class `glmMixture`.

- dispersion:

  The dispersion parameter for the family used. If NULL, it is inferred
  from object.

- ...:

  Additional arguments.

## Value

An object of class `summary.glmMixture` containing:

- call:

  The component from object.

- family:

  The component from object.

- df.residual:

  The residual degrees of freedom.

- coefficients:

  Matrix of coefficients for the outcome model.

- m.coefficients:

  Matrix of coefficients for the mismatch model.

- dispersion:

  Estimated dispersion parameter.

- cov.unscaled:

  The estimated covariance matrix.

- match.prob:

  The posterior match probabilities.

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

summary(fit)
#> 
#> Call:
#> plglm(formula = age_at_death ~ poly(unit_yob, 3, raw = TRUE), 
#>     family = "gaussian", adjustment = adj_object)
#> 
#> Outcome Model Coefficients:
#>                                Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)                     57.7528     0.7554  76.449  < 2e-16 ***
#> poly(unit_yob, 3, raw = TRUE)1 -43.7603     8.6020  -5.087 3.84e-07 ***
#> poly(unit_yob, 3, raw = TRUE)2 114.9040    22.9688   5.003 5.96e-07 ***
#> poly(unit_yob, 3, raw = TRUE)3 -57.1416    16.0604  -3.558 0.000379 ***
#> 
#> Mismatch Model Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)    7.562      3.150   2.400   0.0164 *
#> commf         -6.731      2.839  -2.371   0.0178 *
#> comml         -8.974      3.643  -2.463   0.0138 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 373.1)
#> 
#> Average Correct Match Probability: 0.9506 
#> 
```
