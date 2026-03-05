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
set.seed(124)
linked_df <- data.frame(
  exposure = sample(c("low", "high"), 300, replace = TRUE),
  disease = sample(c("yes", "no"), 300, replace = TRUE)
)

adj <- adjMixture(linked.data = linked_df, m.rate = 0.15)
fit <- plctable(~ exposure + disease, adjustment = adj)

# Extract the variance-covariance matrix of the cell probabilities
vmat <- vcov(fit)
print(vmat)
} # }
```
