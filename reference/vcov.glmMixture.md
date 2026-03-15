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
#> coef (Intercept)                          2.46975902
#> coef poly(unit_yob, 3, raw = TRUE)1     -23.52175723
#> coef poly(unit_yob, 3, raw = TRUE)2      50.32131255
#> coef poly(unit_yob, 3, raw = TRUE)3     -28.84771321
#> dispersion                               13.50918919
#> m.coef (Intercept)                       -0.34449875
#> m.coef commf                              0.09520238
#> m.coef comml                              0.41118201
#>                                     coef poly(unit_yob, 3, raw = TRUE)1
#> coef (Intercept)                                            -23.5217572
#> coef poly(unit_yob, 3, raw = TRUE)1                         330.8987261
#> coef poly(unit_yob, 3, raw = TRUE)2                        -799.9390972
#> coef poly(unit_yob, 3, raw = TRUE)3                         492.3118823
#> dispersion                                                 -176.5971124
#> m.coef (Intercept)                                            4.1178131
#> m.coef commf                                                  0.4020897
#> m.coef comml                                                 -5.6864113
#>                                     coef poly(unit_yob, 3, raw = TRUE)2
#> coef (Intercept)                                              50.321313
#> coef poly(unit_yob, 3, raw = TRUE)1                         -799.939097
#> coef poly(unit_yob, 3, raw = TRUE)2                         2084.357322
#> coef poly(unit_yob, 3, raw = TRUE)3                        -1356.437056
#> dispersion                                                   319.984051
#> m.coef (Intercept)                                            -7.291627
#> m.coef commf                                                  -4.228244
#> m.coef comml                                                  13.515778
#>                                     coef poly(unit_yob, 3, raw = TRUE)3
#> coef (Intercept)                                            -28.8477132
#> coef poly(unit_yob, 3, raw = TRUE)1                         492.3118823
#> coef poly(unit_yob, 3, raw = TRUE)2                       -1356.4370562
#> coef poly(unit_yob, 3, raw = TRUE)3                         931.4044398
#> dispersion                                                    0.1667215
#> m.coef (Intercept)                                            4.5079215
#> m.coef commf                                                  4.3081749
#> m.coef comml                                                -11.3205245
#>                                       dispersion m.coef (Intercept)
#> coef (Intercept)                      13.5091892         -0.3444988
#> coef poly(unit_yob, 3, raw = TRUE)1 -176.5971124          4.1178131
#> coef poly(unit_yob, 3, raw = TRUE)2  319.9840509         -7.2916270
#> coef poly(unit_yob, 3, raw = TRUE)3    0.1667215          4.5079215
#> dispersion                          2018.5407711         15.5338188
#> m.coef (Intercept)                    15.5338188          6.1090919
#> m.coef commf                          -1.3092811         -4.1172117
#> m.coef comml                         -31.3783991         -5.2864219
#>                                     m.coef commf m.coef comml
#> coef (Intercept)                      0.09520238    0.4111820
#> coef poly(unit_yob, 3, raw = TRUE)1   0.40208972   -5.6864113
#> coef poly(unit_yob, 3, raw = TRUE)2  -4.22824390   13.5157781
#> coef poly(unit_yob, 3, raw = TRUE)3   4.30817487  -11.3205245
#> dispersion                           -1.30928107  -31.3783991
#> m.coef (Intercept)                   -4.11721171   -5.2864219
#> m.coef commf                          5.02343397    0.2860827
#> m.coef comml                          0.28608268   10.0709274
```
