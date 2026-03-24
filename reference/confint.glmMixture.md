# Confidence Intervals for glmMixture Objects

Computes Wald confidence intervals for one or more parameters in a
`glmMixture` object.

## Usage

``` r
# S3 method for class 'glmMixture'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `glmMixture`.

- parm:

  A specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  The confidence level required.

- ...:

  Additional arguments (currently ignored).

## Value

A matrix (or vector) with columns giving lower and upper confidence
limits for each parameter.

## Details

The intervals are calculated based on the sandwich variance estimator:
`Estimate +/- z_crit * SE`. For Gaussian and Gamma families, a
t-distribution is used with residual degrees of freedom. For Binomial
and Poisson families, a standard normal distribution is used.

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

confint(fit)
#>                                           2.5 %     97.5 %
#> coef (Intercept)                      54.671426  60.834082
#> coef poly(unit_yob, 3, raw = TRUE)1  -79.426563  -8.093938
#> coef poly(unit_yob, 3, raw = TRUE)2   25.388787 204.419170
#> coef poly(unit_yob, 3, raw = TRUE)3 -116.979978   2.696778
#> dispersion                           285.021138 461.202280
#> m.coef (Intercept)                   -12.407800  -2.715452
#> m.coef commf                           2.336753  11.125784
#> m.coef comml                           2.751873  15.196314
```
