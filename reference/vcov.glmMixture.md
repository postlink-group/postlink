# Extract Variance-Covariance Matrix from a glmMixture Object

Returns the variance-covariance matrix of the main parameters of a
fitted `glmMixture` object. The matrix is estimated using a sandwich
estimator to account for the mixture structure.

## Usage

``` r
# S3 method for class 'glmMixture'
vcov(object, ...)
```

## Arguments

- object:

  An object of class `glmMixture`.

- ...:

  Additional arguments (currently ignored).

## Value

A matrix of the estimated covariances between the parameter estimates.
Row and column names correspond to the parameter names (coefficients,
dispersion, etc.).

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

vcov(fit)
#>                                     coef (Intercept)
#> coef (Intercept)                          2.46975901
#> coef poly(unit_yob, 3, raw = TRUE)1     -23.52175721
#> coef poly(unit_yob, 3, raw = TRUE)2      50.32131252
#> coef poly(unit_yob, 3, raw = TRUE)3     -28.84771319
#> dispersion                               13.50918922
#> m.coef (Intercept)                       -0.34449875
#> m.coef commf                              0.09520238
#> m.coef comml                              0.41118201
#>                                     coef poly(unit_yob, 3, raw = TRUE)1
#> coef (Intercept)                                            -23.5217572
#> coef poly(unit_yob, 3, raw = TRUE)1                         330.8987258
#> coef poly(unit_yob, 3, raw = TRUE)2                        -799.9390964
#> coef poly(unit_yob, 3, raw = TRUE)3                         492.3118816
#> dispersion                                                 -176.5971138
#> m.coef (Intercept)                                            4.1178132
#> m.coef commf                                                  0.4020895
#> m.coef comml                                                 -5.6864112
#>                                     coef poly(unit_yob, 3, raw = TRUE)2
#> coef (Intercept)                                              50.321313
#> coef poly(unit_yob, 3, raw = TRUE)1                         -799.939096
#> coef poly(unit_yob, 3, raw = TRUE)2                         2084.357319
#> coef poly(unit_yob, 3, raw = TRUE)3                        -1356.437054
#> dispersion                                                   319.984054
#> m.coef (Intercept)                                            -7.291627
#> m.coef commf                                                  -4.228243
#> m.coef comml                                                  13.515778
#>                                     coef poly(unit_yob, 3, raw = TRUE)3
#> coef (Intercept)                                            -28.8477132
#> coef poly(unit_yob, 3, raw = TRUE)1                         492.3118816
#> coef poly(unit_yob, 3, raw = TRUE)2                       -1356.4370541
#> coef poly(unit_yob, 3, raw = TRUE)3                         931.4044381
#> dispersion                                                    0.1667189
#> m.coef (Intercept)                                            4.5079214
#> m.coef commf                                                  4.3081748
#> m.coef comml                                                -11.3205243
#>                                       dispersion m.coef (Intercept)
#> coef (Intercept)                      13.5091892         -0.3444988
#> coef poly(unit_yob, 3, raw = TRUE)1 -176.5971138          4.1178132
#> coef poly(unit_yob, 3, raw = TRUE)2  319.9840543         -7.2916273
#> coef poly(unit_yob, 3, raw = TRUE)3    0.1667189          4.5079214
#> dispersion                          2018.5407638         15.5338171
#> m.coef (Intercept)                    15.5338171          6.1090918
#> m.coef commf                          -1.3092796         -4.1172118
#> m.coef comml                         -31.3783985         -5.2864216
#>                                     m.coef commf m.coef comml
#> coef (Intercept)                      0.09520238    0.4111820
#> coef poly(unit_yob, 3, raw = TRUE)1   0.40208952   -5.6864112
#> coef poly(unit_yob, 3, raw = TRUE)2  -4.22824344   13.5157779
#> coef poly(unit_yob, 3, raw = TRUE)3   4.30817477  -11.3205243
#> dispersion                           -1.30927956  -31.3783985
#> m.coef (Intercept)                   -4.11721177   -5.2864216
#> m.coef commf                          5.02343418    0.2860824
#> m.coef comml                          0.28608243   10.0709273
```
