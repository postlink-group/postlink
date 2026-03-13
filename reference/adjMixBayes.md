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
data(lifem)

# lifem data preprocessing
# For computational efficiency in the example, we work with a subset of the lifem data.
lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
lifem_small <- rbind(
  head(subset(lifem, hndlnk == 1), 100),
  head(subset(lifem, hndlnk == 0), 20)
)

# Construct the Bayesian mixture adjustment object
adj_bayes <- adjMixBayes(
  linked.data = lifem_small,
  priors = list(theta = "beta(2, 2)")
)

class(adj_bayes)
#> [1] "adjMixBayes" "adjustment" 
```
