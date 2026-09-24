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
#> coef (Intercept)                           0.5706938
#> coef poly(unit_yob, 3, raw = TRUE)1       -5.2361970
#> coef poly(unit_yob, 3, raw = TRUE)2       11.5372246
#> coef poly(unit_yob, 3, raw = TRUE)3       -7.0462799
#> dispersion                                -0.6942789
#> m.coef (Intercept)                         0.1119627
#> m.coef commf                              -0.0322021
#> m.coef comml                              -0.1508786
#>                                     coef poly(unit_yob, 3, raw = TRUE)1
#> coef (Intercept)                                             -5.2361970
#> coef poly(unit_yob, 3, raw = TRUE)1                          73.9942884
#> coef poly(unit_yob, 3, raw = TRUE)2                        -189.1367509
#> coef poly(unit_yob, 3, raw = TRUE)3                         124.2453630
#> dispersion                                                    9.3269443
#> m.coef (Intercept)                                           -1.1843453
#> m.coef commf                                                 -0.1058184
#> m.coef comml                                                  1.8711390
#>                                     coef poly(unit_yob, 3, raw = TRUE)2
#> coef (Intercept)                                              11.537225
#> coef poly(unit_yob, 3, raw = TRUE)1                         -189.136751
#> coef poly(unit_yob, 3, raw = TRUE)2                          527.564297
#> coef poly(unit_yob, 3, raw = TRUE)3                         -363.681651
#> dispersion                                                   -15.515215
#> m.coef (Intercept)                                             1.233147
#> m.coef commf                                                   1.791217
#> m.coef comml                                                  -3.406404
#>                                     coef poly(unit_yob, 3, raw = TRUE)3
#> coef (Intercept)                                             -7.0462799
#> coef poly(unit_yob, 3, raw = TRUE)1                         124.2453630
#> coef poly(unit_yob, 3, raw = TRUE)2                        -363.6816507
#> coef poly(unit_yob, 3, raw = TRUE)3                         257.9355737
#> dispersion                                                   -0.3502309
#> m.coef (Intercept)                                            0.6959974
#> m.coef commf                                                 -2.5538254
#> m.coef comml                                                  1.1363067
#>                                      dispersion m.coef (Intercept) m.coef commf
#> coef (Intercept)                     -0.6942789          0.1119627   -0.0322021
#> coef poly(unit_yob, 3, raw = TRUE)1   9.3269443         -1.1843453   -0.1058184
#> coef poly(unit_yob, 3, raw = TRUE)2 -15.5152148          1.2331471    1.7912170
#> coef poly(unit_yob, 3, raw = TRUE)3  -0.3502309          0.6959974   -2.5538254
#> dispersion                          120.2448040        -10.7999342    9.7312601
#> m.coef (Intercept)                  -10.7999342          9.9230118   -7.5922150
#> m.coef commf                          9.7312601         -7.5922150    8.0617133
#> m.coef comml                          8.7681912         -8.8896034    3.7141259
#>                                     m.coef comml
#> coef (Intercept)                      -0.1508786
#> coef poly(unit_yob, 3, raw = TRUE)1    1.8711390
#> coef poly(unit_yob, 3, raw = TRUE)2   -3.4064039
#> coef poly(unit_yob, 3, raw = TRUE)3    1.1363067
#> dispersion                             8.7681912
#> m.coef (Intercept)                    -8.8896034
#> m.coef commf                           3.7141259
#> m.coef comml                          13.2713805
```
