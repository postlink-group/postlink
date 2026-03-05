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
set.seed(42)
n <- 100
linked_df <- data.frame(
  x = rnorm(n),
  y = rnorm(n), # In practice, y would contain mismatch errors
  match_score = runif(n, 0.7, 1.0)
)
adj_obj <- adjMixBayes(linked.data = linked_df)

# Implicitly calls print.adjMixBayes()
adj_obj
#> 
#> --- Adjustment Object: Bayesian Mixture (Gutman et al., 2016) ---
#> 
#> * Linked Data:
#>     Observations: 100
#> * Priors:
#>     Status:       None specified. Using symmetric defaults.
#>     Defaults applied during fitting:
#>       GLM:        beta ~ normal(0,5) [binomial: normal(0,2.5)]
#>       Survival:   beta ~ normal(0,5) [weibull: normal(0,2)]
#>       Mix Weight: theta ~ beta(1,1)
#> 
```
