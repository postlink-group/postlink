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
set.seed(125)
n <- 300

# 1. Simulate true categorical data with dependency
exposure <- sample(c("low", "high"), n, replace = TRUE)

# Induce dependency - High exposure -> higher disease probability
prob_disease <- ifelse(exposure == "high", 0.7, 0.3)
true_disease <- ifelse(runif(n) < prob_disease, "yes", "no")

# 2. Induce 15% linkage error
mis_idx <- sample(1:n, size = floor(0.15 * n))
obs_disease <- true_disease
obs_disease[mis_idx] <- sample(obs_disease[mis_idx])

linked_df <- data.frame(exposure = exposure, disease = obs_disease)

# 3. Fit the adjusted contingency table model
adj <- adjMixture(linked.data = linked_df, m.rate = 0.15)
fit <- plctable(~ exposure + disease, adjustment = adj)

# 4. Explicitly call the print method
print(fit)
} # }
```
