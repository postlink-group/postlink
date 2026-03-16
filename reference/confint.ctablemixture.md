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

# 4. Compute confidence intervals
# 95% CI for all cell probabilities
confint(fit)

# 90% CI for specific cells by name
confint(fit, parm = c("(low, yes)", "(high, no)"), level = 0.90)
} # }
```
