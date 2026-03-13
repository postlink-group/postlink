# Print pooled regression results

Print pooled regression results

## Usage

``` r
# S3 method for class 'mi_link_pool_glm'
print(x, digits = max(3L, getOption("digits") - 2L), ...)
```

## Arguments

- x:

  An object of class `mi_link_pool_glm`, typically returned by
  [`mi_with()`](https://postlink-group.github.io/postlink/reference/mi_with.md)
  for a `glmMixBayes` fit.

- digits:

  the number of significant digits to print.

- ...:

  further arguments (unused).

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

pooled_fit <- mi_with(
  object = fit,
  data = lifem_small,
  formula = age_at_death ~ poly(unit_yob, 3, raw = TRUE),
  family = gaussian()
)

print(pooled_fit, digits = 4)
} # }
```
