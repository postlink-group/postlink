# Print a glmMixBayes model object

Print a glmMixBayes model object

## Usage

``` r
# S3 method for class 'glmMixBayes'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `glmMixBayes`.

- digits:

  Minimum number of significant digits to show.

- ...:

  Further arguments (unused).

## Value

The input `x`, invisibly.

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

print(fit)
} # }
```
