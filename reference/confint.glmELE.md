# Confidence Intervals for `glmELE` Objects

Computes Wald confidence intervals for one or more parameters in a
`glmELE` object.

## Usage

``` r
# S3 method for class 'glmELE'
confint(object, parm, level = 0.95, weight.matrix = NULL, ...)
```

## Arguments

- object:

  An object of class `"glmELE"`.

- parm:

  A specification of which parameters are to be given confidence
  intervals, either a vector of numbers or names. If missing, all
  parameters are considered.

- level:

  The confidence level required.

- weight.matrix:

  Character string specifying the weighting method to use. Defaults to
  the first method found.

- ...:

  Additional arguments passed to methods.

## Value

A matrix (or vector) with columns giving lower and upper confidence
limits for each parameter.

## Examples

``` r
data(brfss, package = "postlink")

adj_object <- adjELE(linked.data = brfss,
                     m.rate = unique(brfss$m.rate),
                     blocks = imonth,
                     weight.matrix = "BLUE")

fit <- plglm(Weight ~ Height + Physhlth + Menthlth + Exerany,
             family = "gaussian", adjustment = adj_object)

confint(fit)
#>                     2.5 %      97.5 %
#> (Intercept) -1.594774e+02 -78.3262026
#> Height       5.121013e-01   0.6639389
#> Physhlth     9.808049e-03   0.6347342
#> Menthlth    -2.309497e-01   0.4939633
#> Exerany     -8.562848e+00  -0.3385568
```
