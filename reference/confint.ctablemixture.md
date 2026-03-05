# Confidence Intervals for Adjusted Cell Probabilities

Computes Wald-type confidence intervals for the estimated cell
probabilities of the correctly matched population.

## Usage

``` r
# S3 method for class 'ctableMixture'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `ctableMixture`.

- parm:

  A specification of which parameters are to be given confidence
  intervals. If missing, all parameters (cells) are considered.

- level:

  The confidence level required. Defaults to 0.95.

- ...:

  Additional arguments (currently ignored).

## Value

A matrix with columns giving lower and upper confidence limits for each
parameter.

## Details

The intervals are calculated using the standard error estimates derived
from
[`vcov.ctableMixture`](https://postlink-group.github.io/postlink/reference/vcov.ctablemixture.md).
The lower and upper bounds are truncated at 0 and 1, respectively, to
ensure valid probability estimates.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(125)
linked_df <- data.frame(
  exposure = sample(c("low", "high"), 300, replace = TRUE),
  disease = sample(c("yes", "no"), 300, replace = TRUE)
)

adj <- adjMixture(linked.data = linked_df, m.rate = 0.15)
fit <- plctable(~ exposure + disease, adjustment = adj)

# Compute 95% confidence intervals for all cell probabilities
confint(fit)

# Compute 90% confidence intervals for specific cells by name
confint(fit, parm = c("(low, yes)", "(high, no)"), level = 0.90)
} # }
```
