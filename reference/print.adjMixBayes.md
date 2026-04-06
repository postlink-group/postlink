# Print Method for `adjMixBayes` Objects

Provides a concise summary of the Bayesian adjustment object created by
[`adjMixBayes`](https://postlink-group.github.io/postlink/reference/adjMixBayes.md).

## Usage

``` r
# S3 method for class 'adjMixBayes'
print(x, ...)
```

## Arguments

- x:

  An object of class `adjMixBayes`.

- ...:

  Additional arguments passed to methods.

## Value

Invisibly returns `x`.

## Details

This method inspects the reference-based data environment to report the
number of linked records without copying the full dataset. It safely
handles cases where the linked data is unspecified (NULL). It also
prints the user-specified priors or outlines the defaults that will be
used.

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

adj_obj <- adjMixBayes(
  linked.data = lifem_small,
  priors = list(theta = "beta(2, 2)")
)

# Implicitly calls print.adjMixBayes()
adj_obj
#> 
#> --- Adjustment Object: Bayesian Mixture (Gutman et al., 2016) ---
#> 
#> * Linked Data:
#>     Observations: 120
#> * Priors:
#>     User-specified overrides:
#>       theta      : beta(2, 2)
#>     (Unspecified parameters will use defaults below)
#>     Defaults applied during fitting:
#>       Intercept:  intercept ~ normal(0,10)
#>       GLM Slopes: beta ~ normal(0,5) [binomial: normal(0,2.5)]
#>       Surv Slopes: beta ~ normal(0,5) [weibull: normal(0,2)]
#>       Mix Weight:  theta ~ beta(1,1)
#> 
```
