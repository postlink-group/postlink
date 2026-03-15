# Print pooled regression results

Print pooled regression results

## Usage

``` r
# S3 method for class 'mi_link_pool_glm'
print(x, digits = max(3L, getOption("digits") - 2L), ...)
```

## Arguments

- x:

  An object of class `mi_link_pool_glm`, typically returned by
  [`mi_with()`](https://postlink-group.github.io/postlink/reference/mi_with.md)
  for a `glmMixBayes` fit.

- digits:

  the number of significant digits to print.

- ...:

  further arguments (unused).

## Value

The input `x`, invisibly.

## Examples

``` r
# \donttest{
# Simulate data
set.seed(607)
n <- 100
linked_data <- data.frame(x1 = rnorm(n), y = rnorm(n))
X <- model.matrix(~ x1, data = linked_data)

# Fit the mixture model
fit <- glmMixBayes(
  X = X, y = linked_data$y, family = "gaussian",
  control = list(iterations = 150, burnin.iterations = 50, seed = 607)
)
#> 
#> SAMPLING FOR MODEL 'glmMixBayes_gaussian' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.58 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 7
#> Chain 1:            adapt_window = 38
#> Chain 1:            term_buffer = 5
#> Chain 1: 
#> Chain 1: Iteration:   1 / 150 [  0%]  (Warmup)
#> Chain 1: Iteration:  15 / 150 [ 10%]  (Warmup)
#> Chain 1: Iteration:  30 / 150 [ 20%]  (Warmup)
#> Chain 1: Iteration:  45 / 150 [ 30%]  (Warmup)
#> Chain 1: Iteration:  51 / 150 [ 34%]  (Sampling)
#> Chain 1: Iteration:  65 / 150 [ 43%]  (Sampling)
#> Chain 1: Iteration:  80 / 150 [ 53%]  (Sampling)
#> Chain 1: Iteration:  95 / 150 [ 63%]  (Sampling)
#> Chain 1: Iteration: 110 / 150 [ 73%]  (Sampling)
#> Chain 1: Iteration: 125 / 150 [ 83%]  (Sampling)
#> Chain 1: Iteration: 140 / 150 [ 93%]  (Sampling)
#> Chain 1: Iteration: 150 / 150 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.14 seconds (Warm-up)
#> Chain 1:                0.117 seconds (Sampling)
#> Chain 1:                0.257 seconds (Total)
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
#> Global label swap performed: label 2 dominates label 1.
#> 
#>     ......................................................................................
#>     . Method                         Time (sec)           Status                         . 
#>     ......................................................................................
#>     . ECR-ITERATIVE-1                0.112                Converged (3 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

# Pool regression fits across posterior draws
pooled_fit <- mi_with(
  object = fit,
  data = linked_data,
  formula = y ~ x1,
  family = gaussian()
)

# Explicitly test the print method
print(pooled_fit, digits = 4)
#> Pooled regression results across posterior match classifications:
#>   Retained imputations (m): 100 
#> 
#>             Estimate Std.Error  CI.lwr CI.upr       df
#> (Intercept)   0.1001    0.1667 -0.2281 0.4282 302.8298
#> x1            0.0138    0.1447 -0.2702 0.2979 630.7583
# }
```
