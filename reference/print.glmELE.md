# Print a `glmELE` Object

Prints the function call and the estimated coefficient matrices from a
fitted `glmELE` object.

## Usage

``` r
# S3 method for class 'glmELE'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"glmELE"`.

- digits:

  The number of significant digits to print. Defaults to
  `max(3L, getOption("digits") - 3L)`.

- ...:

  Additional arguments passed to methods.

## Value

Invisibly returns the input object `x`.

## Examples

``` r
data(brfss, package = "postlink")

adj_object <- adjELE(linked.data = brfss,
                     m.rate = unique(brfss$m.rate),
                     blocks = imonth,
                     weight.matrix = "BLUE")

fit <- plglm(Weight ~ Height + Physhlth + Menthlth + Exerany,
             family = "gaussian", adjustment = adj_object)

print(fit)
#> 
#> Call:
#> plglm(formula = Weight ~ Height + Physhlth + Menthlth + Exerany, 
#>     family = "gaussian", adjustment = adj_object)
#> 
#> Coefficients:
#>       (Intercept)  Height     Physhlth   Menthlth   Exerany  
#> BLUE  -118.9018       0.5880     0.3223     0.1315    -4.4507
#> 
```
