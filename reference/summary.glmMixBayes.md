# Summary method for glmMixBayes models

Summary method for glmMixBayes models

## Usage

``` r
# S3 method for class 'glmMixBayes'
summary(object, ...)
```

## Arguments

- object:

  An object of class `glmMixBayes`.

- ...:

  Not used.

## Value

An object of class `"summary.glmMixBayes"`, which is printed with a
custom method.

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
#> Chain 1: Gradient evaluation took 4.9e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.49 seconds.
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
#> Chain 1:  Elapsed Time: 0.599 seconds (Warm-up)
#> Chain 1:                0.771 seconds (Sampling)
#> Chain 1:                1.37 seconds (Total)
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
#>     . ECR-ITERATIVE-1                0.076                Converged (2 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

summary(fit)
#> Call:
#> plglm(formula = age_at_death ~ poly(unit_yob, 3, raw = TRUE), 
#>     family = "gaussian", adjustment = adj, control = list(iterations = 200, 
#>         burnin.iterations = 100, seed = 123))
#>  
#> Family:
#> gaussian
#>  
#> (Component 1 = Correct-match):
#> Outcome Model Coefficients:
#>                                Estimates Std. Error  2.5 % 97.5 %
#> (Intercept)                       44.890      2.095 40.870  49.12
#> poly(unit_yob, 3, raw = TRUE)1    11.604      4.324  3.284  18.83
#> poly(unit_yob, 3, raw = TRUE)2     8.059      4.309 -2.128  15.62
#> poly(unit_yob, 3, raw = TRUE)3     8.326      5.146 -3.542  18.75
#>  
#> Dispersion:
#>   Estimate  Std. Error
#>   350.8      64.9     
#> 
#> (Component 2 = Incorrect-match):
#> Outcome Model Coefficients:
#>                                Estimates Std. Error   2.5 % 97.5 %
#> (Intercept)                       2.5058     4.6753 -6.1045 10.903
#> poly(unit_yob, 3, raw = TRUE)1    1.4558     5.1252 -7.6492 10.525
#> poly(unit_yob, 3, raw = TRUE)2    0.9269     4.6265 -7.9656  8.512
#> poly(unit_yob, 3, raw = TRUE)3    1.2907     4.6886 -7.4155  9.298
#>  
#> Dispersion:
#>   Estimate  Std. Error
#>   130.6     373.2     
#> 
# }
```
