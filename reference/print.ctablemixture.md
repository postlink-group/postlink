# Print Method for Adjusted Contingency Tables

Prints the estimated contingency table (corrected for linkage error) and
a summary of the adjustment parameters used by the mixture model.

## Usage

``` r
# S3 method for class 'ctableMixture'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `ctableMixture`.

- digits:

  Integer; the number of significant digits to use when printing numeric
  values. Defaults to 3.

- ...:

  Additional arguments passed to
  [`print.default`](https://rdrr.io/r/base/print.default.html).

## Value

The argument `x`, invisibly.

## See also

[`plctable`](https://postlink-group.github.io/postlink/reference/plctable.md),
[`ctableMixture`](https://postlink-group.github.io/postlink/reference/ctableMixture.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fast simulation of linked data
set.seed(123)
linked_df <- data.frame(
  exposure = sample(c("low", "high"), 300, replace = TRUE),
  disease = sample(c("yes", "no"), 300, replace = TRUE)
)

# Specify adjustment and fit the model
adj <- adjMixture(linked.data = linked_df, m.rate = 0.15)
fit <- plctable(~ exposure + disease, adjustment = adj)

# Explicitly call the print method
print(fit)
} # }
```
