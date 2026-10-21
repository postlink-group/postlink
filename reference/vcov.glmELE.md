# Extract Variance-Covariance Matrix from a `glmELE` Object

Extracts the variance-covariance matrix of the main parameters for a
specific weighting method.

## Usage

``` r
# S3 method for class 'glmELE'
vcov(object, weight.matrix = NULL, ...)
```

## Arguments

- object:

  An object of class `"glmELE"`.

- weight.matrix:

  Character string specifying which weighting method to return. Defaults
  to the first method found in the object.

- ...:

  Additional arguments passed to methods.

## Value

A matrix of the estimated covariances between the parameter estimates.

## Examples

``` r
data(brfss, package = "postlink")

adj_object <- adjELE(linked.data = brfss,
                     m.rate = unique(brfss$m.rate),
                     blocks = imonth,
                     weight.matrix = "BLUE")

fit <- plglm(Weight ~ Height + Physhlth + Menthlth + Exerany,
             family = "gaussian", adjustment = adj_object)

vcov(fit)
#>              (Intercept)        Height      Physhlth      Menthlth
#> (Intercept) 428.04649877 -7.700850e-01 -1.933555e-01 -8.019826e-02
#> Height       -0.77008499  1.498509e-03  3.004707e-05  6.287311e-05
#> Physhlth     -0.19335545  3.004707e-05  2.538390e-02 -1.142878e-02
#> Menthlth     -0.08019826  6.287311e-05 -1.142878e-02  3.415645e-02
#> Exerany     -10.82311640 -1.651295e-03  4.033932e-02 -2.424824e-03
#>                   Exerany
#> (Intercept) -10.823116402
#> Height       -0.001651295
#> Physhlth      0.040339315
#> Menthlth     -0.002424824
#> Exerany       4.396407143
```
