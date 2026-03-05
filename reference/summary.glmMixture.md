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
#> (Intercept)                      57.753      1.572  36.749   <2e-16 ***
#> poly(unit_yob, 3, raw = TRUE)1  -43.760     18.191  -2.406   0.0162 *  
#> poly(unit_yob, 3, raw = TRUE)2  114.904     45.655   2.517   0.0119 *  
#> poly(unit_yob, 3, raw = TRUE)3  -57.142     30.519  -1.872   0.0613 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Mismatch Model Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept)   -7.562      2.472  -3.059  0.00222 **
#> commf          6.731      2.241   3.003  0.00267 **
#> comml          8.974      3.173   2.828  0.00469 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 373.1)
#> 
#> Average Correct Match Probability: 0.951 
#> 
```
