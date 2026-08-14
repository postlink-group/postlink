# Unsupported Generics for postlink Models

Intercepts standard likelihood and residual-based generics that are not
directly applicable to post-linkage models because the true match status
remains latent.

## Usage

``` r
# S3 method for class 'plmodel'
logLik(object, ...)

# S3 method for class 'plmodel'
profile(fitted, ...)

# S3 method for class 'plmodel'
anova(object, ...)

# S3 method for class 'plmodel'
extractAIC(fit, scale, k = 2, ...)

# S3 method for class 'plglm'
cooks.distance(model, ...)

# S3 method for class 'plglm'
rstudent(model, ...)
```

## Arguments

- object:

  A fitted model object.

- ...:

  Additional arguments ignored by these methods.

- fitted:

  A fitted model object.

- fit:

  A fitted model object.

- scale:

  Optional numeric specifying the scale parameter.

- k:

  Numeric specifying the "weight" of the equivalent degrees of freedom.

- model:

  A fitted model object.
