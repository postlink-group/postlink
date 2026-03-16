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

# 4. Generate the detailed summary object
sum_fit <- summary(fit)

# 5. Access specific components of the summary
print(sum_fit$coefficients)
print(sum_fit$chisq)
} # }
```
