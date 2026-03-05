# Contingency Table Analysis with Mixture-Based Adjustment

Estimates the cell probabilities of a two-way contingency table in the
presence of linkage errors. The function implements the methodology
described in Slawski et al. (2025), modeling the observed table as a
mixture of correctly matched records (following a saturated model) and
mismatched records (assumed to follow an independence model).

## Usage

``` r
ctableMixture(tab, m.rate, control = list(), ...)
```

## Arguments

- tab:

  A numeric matrix or table of counts representing the observed two-way
  contingency table.

- m.rate:

  A numeric value between 0 and 1 indicating the assumed rate of
  mismatched records in the data.

- control:

  A list of control parameters. See `...` for details.

- ...:

  Additional control arguments. If not provided in `control`, these will
  be used. Supported arguments:

  - `max.iter`: Integer. Maximum number of EM iterations (default:
    1000).

  - `tol`: Numeric. Convergence tolerance for the negative
    log-likelihood (default: 1e-6).

## Value

A list of results:

- `phat`: Matrix of estimated cell probabilities for the correctly
  matched population (the target parameter).

- `phat0`: Matrix of estimated cell probabilities for the mismatched
  population (independence model).

- `var`: Estimated variance-covariance matrix of the estimators in
  `phat`.

- `ftable`: The estimated contingency table of counts for the correctly
  matched population (adjusted for bias).

- `objective`: The final value of the negative log-likelihood.

- `converged`: Logical indicating if the algorithm converged within
  `max.iter`.

## Details

In the absence of linkage errors, the standard estimator for cell
probabilities is the matrix of observed relative frequencies. When
linkage errors are present (specifically mismatches), this estimator is
biased.

This function corrects for this bias using an Expectation-Maximization
(EM) algorithm. The observed data is modeled as a mixture: \$\$P\_{obs}
= (1 - \alpha) P\_{sat} + \alpha P\_{ind}\$\$ where:

- \\\alpha\\ (`m.rate`) is the mismatch rate (fixed/known).

- \\P\_{sat}\\ is the distribution of correct matches (saturated model).

- \\P\_{ind}\\ is the distribution of mismatches (independence model).

The algorithm iteratively updates the posterior probability of a record
being a correct match (E-step) and the estimates for the saturated and
independence distributions (M-step).

## References

Slawski, M., West, B. T., Bukke, P., Wang, Z., Diao, G., & Ben-David, E.
(2025). A general framework for regression with mismatched data based on
mixture modelling. *Journal of the Royal Statistical Society Series A*.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate Synthetic Data
set.seed(1234)
K <- 3; L <- 4
n <- 1000
# Define true probabilities for a KxL table
cellprobs <- c(0.18, 0.05, 0.03, 0.04, 0.02, 0.14,
               0.02, 0.02, 0.10, 0.21, 0.15, 0.04)
matrix_probs <- matrix(cellprobs, nrow = K, ncol = L)

# Generate multinomial counts
dat <- stats::rmultinom(n = n, size = 1, prob = cellprobs)
obs_idx <- apply(dat, 2, function(x) which(x == 1))
X <- ceiling(obs_idx / L) # Row indices
Y <- (obs_idx %% L); Y[Y == 0] <- L # Col indices

# Introduce Linkage Error (Mismatches)
alpha <- 0.20 # 20% mismatch rate
n_mismatch <- round(n * alpha)
Y_perm <- Y
# Shuffle the first n_mismatch Y values to break dependence
Y_perm[1:n_mismatch] <- sample(Y[1:n_mismatch])

# Create Observed Table (with error)
tab_obs <- table(X, Y_perm)

# Apply Adjustment Method
fit <- ctableMixture(tab = tab_obs, m.rate = alpha)

# Inspect Results
print(fit$converged)
# Compare estimated Correct Counts vs True Counts (approx)
print(round(fit$ftable))
print(round(table(X, Y))) # True table without errors
} # }
```
