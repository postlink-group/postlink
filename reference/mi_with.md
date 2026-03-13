# Pool parameter estimates across posterior draws

Generic function for pooling parameter estimates from Bayesian mixture
models using posterior draws of the latent match classifications.

## Usage

``` r
mi_with(object, ...)
```

## Arguments

- object:

  A fitted Bayesian mixture model object.

- ...:

  Additional arguments passed to methods.

## Value

A pooled model object combining parameter estimates across posterior
component allocations.

## Examples

``` r
if (FALSE) { # \dontrun{
# mi_with() is a generic function for posterior allocation–based pooling.
# See ?mi_with.glmMixBayes for a complete example illustrating its use
# with Bayesian GLM mixture models.
} # }
```
