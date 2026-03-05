# Summarize a `glmELE` Object

Summarizes the results from a `glmELE` fit, providing coefficient
estimates, standard errors, test statistics, and p-values for each
weighting method used.

## Usage

``` r
# S3 method for class 'glmELE'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"glmELE"`.

- ...:

  Additional arguments passed to methods.

## Value

An object of class `"summary.glmELE"`, which is a list containing:

- call:

  The matched call.

- family:

  The family object used.

- coefficients:

  A list of matrices, one per weighting method, containing estimates,
  SEs, t/z values, and p-values.

- dispersion:

  The estimated dispersion parameter(s).

- deviance:

  The deviance of the fitted model.

- df.residual:

  The residual degrees of freedom.

## Examples

``` r
data(brfss, package = "postlink")

adj_object <- adjELE(linked.data = brfss,
                     m.rate = unique(brfss$m.rate),
                     blocks = imonth,
                     weight.matrix = "BLUE")

fit <- plglm(Weight ~ Height + Physhlth + Menthlth + Exerany,
             family = "gaussian", adjustment = adj_object)

summary(fit)
#> 
#> Call:
#> plglm(formula = Weight ~ Height + Physhlth + Menthlth + Exerany, 
#>     family = "gaussian", adjustment = adj_object)
#> 
#> Family: gaussian
#> 
#> --- Weighting Method: BLUE ---
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -118.90181   20.68928  -5.747 1.05e-08 ***
#> Height         0.58802    0.03871  15.190  < 2e-16 ***
#> Physhlth       0.32227    0.15932   2.023   0.0432 *  
#> Menthlth       0.13151    0.18481   0.712   0.4768    
#> Exerany       -4.45070    2.09676  -2.123   0.0339 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Dispersion parameter for gaussian family: 1434
#> 
```
