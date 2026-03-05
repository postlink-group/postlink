# Summary Method for Adjusted Contingency Tables

Provides a detailed summary of the `ctableMixture` model fit, including
the estimated cell probabilities with standard errors, convergence
status, and a Chi-squared test of independence performed on the adjusted
counts.

## Usage

``` r
# S3 method for class 'ctableMixture'
summary(object, ...)
```

## Arguments

- object:

  An object of class `ctableMixture`.

- ...:

  Additional arguments (currently ignored).

## Value

An object of class `summary.ctableMixture` containing:

- call:

  The function call.

- m.rate:

  The assumed mismatch rate.

- ftable:

  The estimated contingency table of correctly matched counts.

- coefficients:

  A matrix containing estimates, standard errors, z-values, and p-values
  for cell probabilities.

- chisq:

  The result of a Pearson's Chi-squared test on the adjusted table.

- converged:

  Logical indicating if the EM algorithm converged.

- iterations:

  Number of iterations performed.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(126)
linked_df <- data.frame(
  exposure = sample(c("low", "high"), 300, replace = TRUE),
  disease = sample(c("yes", "no"), 300, replace = TRUE)
)

adj <- adjMixture(linked.data = linked_df, m.rate = 0.15)
fit <- plctable(~ exposure + disease, adjustment = adj)

# Generate the detailed summary object
sum_fit <- summary(fit)

# Access specific components of the summary
print(sum_fit$coefficients)
print(sum_fit$chisq)
} # }
```
