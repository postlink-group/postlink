# Extract Variance-Covariance Matrix from ctableMixture Objects

Extracts the estimated variance-covariance matrix of the cell
probabilities from a fitted `ctableMixture` object. The variance is
estimated using the observed information matrix (via the Hessian of the
mixture log-likelihood).

## Usage

``` r
# S3 method for class 'ctableMixture'
vcov(object, ...)
```

## Arguments

- object:

  An object of class `ctableMixture`.

- ...:

  Additional arguments (currently ignored).

## Value

A matrix of the estimated covariances between the cell probability
estimates. The row and column names correspond to the cells of the table
in row-major order (e.g., "(Row1, Col1)", "(Row1, Col2)", ...).

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

# 4. Extract the variance-covariance matrix of the cell probabilities
vmat <- vcov(fit)
print(vmat)
} # }
```
