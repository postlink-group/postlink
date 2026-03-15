# Summary method for glmMixBayes models

Summary method for glmMixBayes models

## Usage

``` r
# S3 method for class 'glmMixBayes'
summary(object, ...)
```

## Arguments

- object:

  An object of class `glmMixBayes`.

- ...:

  Not used.

## Value

An object of class `"summary.glmMixBayes"`, which is printed with a
custom method.

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

adj <- adjMixBayes(
  linked.data = lifem_small,
  priors = list(theta = "beta(2, 2)")
)

fit <- plglm(
  age_at_death ~ poly(unit_yob, 3, raw = TRUE),
  family = "gaussian",
  adjustment = adj,
  control = list(
    iterations = 200,
    burnin.iterations = 100,
    seed = 123
  )
)

summary(fit)
} # }
```
