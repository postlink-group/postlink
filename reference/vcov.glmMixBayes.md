# Posterior covariance matrix for glmMixBayes coefficients

Posterior covariance matrix for glmMixBayes coefficients

## Usage

``` r
# S3 method for class 'glmMixBayes'
vcov(object, ...)
```

## Arguments

- object:

  A `glmMixBayes` model object.

- ...:

  Not used.

## Value

Posterior covariance matrix of the regression coefficients for the
primary (correct-match) component.

## Examples

``` r
if (FALSE) { # \dontrun{
data(lifem)

# lifem data preprocessing
# For computational efficiency in the example, we work with a subset of the lifem data.
lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
lifem_small <- rbind(
  head(subset(lifem, hndlnk == 1), 100),
  head(subset(lifem, hndlnk == 0), 20)
)

x <- cbind(1, poly(lifem_small$unit_yob, 3, raw = TRUE))
y <- lifem_small$age_at_death

fit <- glmMixBayes(
  X = x,
  y = y,
  family = "gaussian",
  control = list(
    iterations = 200,
    burnin.iterations = 100,
    seed = 123
  )
)

vcov(fit)
} # }
```
