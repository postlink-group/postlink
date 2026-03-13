# Pooling regression fits across posterior draws of correct-match classifications

Use posterior draws of the latent match indicators from
[`glmMixBayes()`](https://postlink-group.github.io/postlink/reference/glmMixBayes.md)
to repeatedly identify which records are treated as correct matches,
refit the requested regression model on those records, and pool the
resulting estimates.

Each retained posterior draw defines one subset of records classified as
correct matches. The function fits the specified
[`lm()`](https://rdrr.io/r/stats/lm.html) or
[`glm()`](https://rdrr.io/r/stats/glm.html) model to that subset,
extracts the estimated coefficients and their covariance matrix, and
combines the results across draws using multiple-imputation pooling
rules.

## Usage

``` r
# S3 method for class 'glmMixBayes'
mi_with(
  object,
  data,
  formula,
  family = NULL,
  min_n = NULL,
  quietly = TRUE,
  ...
)
```

## Arguments

- object:

  A `glmMixBayes` model object containing posterior draws of the latent
  match indicators.

- data:

  A data.frame with all candidate records in the same row order as used
  in the model.

- formula:

  Model formula for refitting on each draw (required).

- family:

  A [`stats::family()`](https://rdrr.io/r/stats/family.html) object for
  the refitted model. If not supplied, the function chooses a default
  family based on `object$family`.

- min_n:

  Minimum number of records required to fit the model for a given
  posterior draw. The default is `p + 1`, where `p` is the number of
  columns in the model matrix.

- quietly:

  If `TRUE`, draws that lead to fitting errors are skipped without
  printing the full error message.

- ...:

  Additional arguments passed through (currently unused).

## Value

An object of class `c("mi_link_pool_glm", "mi_link_pool")` containing
pooled coefficient estimates, standard errors, confidence intervals, and
related summary information.

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

print(pooled_fit)
} # }
```
