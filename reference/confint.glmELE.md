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
#> (Intercept) -1.591870e+02 -78.6166068
#> Height       5.126909e-01   0.6633492
#> Physhlth     9.775369e-03   0.6347669
#> Menthlth    -2.309039e-01   0.4939175
#> Exerany     -8.566838e+00  -0.3345667
```
