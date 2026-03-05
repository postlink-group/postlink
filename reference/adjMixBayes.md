# Secondary Analysis Constructor Based on Bayesian Mixture Modeling

Specifies the linked data and information on the underlying record
linkage process for a mismatch error adjustment using a Bayesian
framework based on mixture modeling as developed by Gutman et al.
(2016). This framework uses a mixture model for pairs of linked records
whose two components reflect distributions conditional on match status,
i.e., correct match or mismatch. Posterior inference is carried out via
data augmentation or multiple imputation.

## Usage

``` r
adjMixBayes(linked.data = NULL, priors = NULL)
```

## Arguments

- linked.data:

  A data.frame containing the linked dataset.

- priors:

  A named `list` (or `NULL`) of prior specifications. Because the Stan
  models are pre-compiled, these strings are parsed into numeric
  hyperparameters and passed to the model's data block. Any missing
  entries are automatically filled with symmetric defaults dynamically
  during the model fitting phase.

## Value

An object of class `c("adjMixBayes", "adjustment")`. To minimize memory
overhead, the underlying `linked.data` is stored by reference within an
environment inside this object.

## Details

Explicit provision of `linked.data` is strongly recommended for
reproducibility and to ensure the adjustment object fully encapsulates
the necessary data for downstream model fitting.

## References

Gutman, R., Sammartino, C., Green, T., & Montague, B. (2016). Error
adjustments for file linking methods using encrypted unique client
identifier (eUCI) with application to recently released prisoners who
are HIV+. *Statistics in Medicine*, 35(1), 115–129.
[doi:10.1002/sim.6586](https://doi.org/10.1002/sim.6586)

## See also

- [`plglm()`](https://postlink-group.github.io/postlink/reference/plglm.md)
  for generalized linear regression modeling

- [`plsurvreg()`](https://postlink-group.github.io/postlink/reference/plsurvreg.md)
  for parametric survival modeling

## Examples

``` r
# Simulate true data
set.seed(123)
n <- 200
x <- rnorm(n)
true_y <- 1.5 + 2 * x + rnorm(n)

# Induce linkage mismatch errors
# Assume we have a match probability for each record
match_prob <- rbeta(n, 8, 2)
is_mismatch <- rbinom(n, 1, 1 - match_prob)

obs_y <- true_y
mismatch_idx <- which(is_mismatch == 1)
if(length(mismatch_idx) > 1) {
  obs_y[mismatch_idx] <- sample(obs_y[mismatch_idx])
}
linked_data <- data.frame(y = obs_y, x = x, match_prob = match_prob)

# Construct the Bayesian mixture adjustment object
adj_bayes <- adjMixBayes(
  linked.data = linked_data,
  priors = list(theta = "beta(2, 2)") # Optional: Override default
)

class(adj_bayes)
#> [1] "adjMixBayes" "adjustment" 
```
