# Posterior covariance matrix for glmMixBayes coefficients

Posterior covariance matrix for glmMixBayes coefficients

## Usage

``` r
# S3 method for class 'glmMixBayes'
vcov(object, ...)
```

## Arguments

- object:

  A `glmMixBayes` model object.

- ...:

  Not used.

## Value

Posterior covariance matrix of the regression coefficients for component
1 (the correct-match component).

## Examples

``` r
# \donttest{
data(lifem)

# lifem data preprocessing
# For computational efficiency in the example, we work with a subset of the lifem data.
lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
lifem_small <- rbind(
  head(subset(lifem, hndlnk == 1), 100),
  head(subset(lifem, hndlnk == 0), 20)
)

x <- cbind(1, poly(lifem_small$unit_yob, 3, raw = TRUE))
y <- lifem_small$age_at_death

adj <- adjMixBayes(
  linked.data = lifem_small,
  priors = list(theta = "beta(2, 2)")
)

fit <- plglm(
  age_at_death ~ poly(unit_yob, 3, raw = TRUE),
  family = "gaussian",
  adjustment = adj,
  control = list(
    iterations = 200,
    burnin.iterations = 100,
    seed = 123
  )
)
#> 
#> SAMPLING FOR MODEL 'glmMixBayes_gaussian' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.47 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 15
#> Chain 1:            adapt_window = 75
#> Chain 1:            term_buffer = 10
#> Chain 1: 
#> Chain 1: Iteration:   1 / 200 [  0%]  (Warmup)
#> Chain 1: Iteration:  20 / 200 [ 10%]  (Warmup)
#> Chain 1: Iteration:  40 / 200 [ 20%]  (Warmup)
#> Chain 1: Iteration:  60 / 200 [ 30%]  (Warmup)
#> Chain 1: Iteration:  80 / 200 [ 40%]  (Warmup)
#> Chain 1: Iteration: 100 / 200 [ 50%]  (Warmup)
#> Chain 1: Iteration: 101 / 200 [ 50%]  (Sampling)
#> Chain 1: Iteration: 120 / 200 [ 60%]  (Sampling)
#> Chain 1: Iteration: 140 / 200 [ 70%]  (Sampling)
#> Chain 1: Iteration: 160 / 200 [ 80%]  (Sampling)
#> Chain 1: Iteration: 180 / 200 [ 90%]  (Sampling)
#> Chain 1: Iteration: 200 / 200 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.64 seconds (Warm-up)
#> Chain 1:                0.762 seconds (Sampling)
#> Chain 1:                1.402 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> 
#>     ......................................................................................
#>     . Method                         Time (sec)           Status                         . 
#>     ......................................................................................
#>     . ECR-ITERATIVE-1                0.079                Converged (2 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

vcov(fit)
#>                                (Intercept) poly(unit_yob, 3, raw = TRUE)1
#> (Intercept)                       4.646952                      -3.037748
#> poly(unit_yob, 3, raw = TRUE)1   -3.037748                      18.637965
#> poly(unit_yob, 3, raw = TRUE)2   -2.332704                      -6.143134
#> poly(unit_yob, 3, raw = TRUE)3    1.676478                      -6.518406
#>                                poly(unit_yob, 3, raw = TRUE)2
#> (Intercept)                                         -2.332704
#> poly(unit_yob, 3, raw = TRUE)1                      -6.143134
#> poly(unit_yob, 3, raw = TRUE)2                      21.137074
#> poly(unit_yob, 3, raw = TRUE)3                      -8.875787
#>                                poly(unit_yob, 3, raw = TRUE)3
#> (Intercept)                                          1.676478
#> poly(unit_yob, 3, raw = TRUE)1                      -6.518406
#> poly(unit_yob, 3, raw = TRUE)2                      -8.875787
#> poly(unit_yob, 3, raw = TRUE)3                      20.120496
# }
```
