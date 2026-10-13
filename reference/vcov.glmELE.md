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
#> (Intercept) 421.94126975 -7.581319e-01 -1.970479e-01 -4.875877e-02
#> Height       -0.75813191  1.475323e-03  3.599829e-05  7.049590e-06
#> Physhlth     -0.19704793  3.599829e-05  2.538921e-02 -1.144378e-02
#> Menthlth     -0.04875877  7.049590e-06 -1.144378e-02  3.414782e-02
#> Exerany     -10.86577018 -1.611498e-03  4.057175e-02 -3.363071e-03
#>                   Exerany
#> (Intercept) -10.865770184
#> Height       -0.001611498
#> Physhlth      0.040571746
#> Menthlth     -0.003363071
#> Exerany       4.404943103
```
