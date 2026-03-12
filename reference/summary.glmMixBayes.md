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
# Simulate data
set.seed(602)
n <- 100
X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("Intercept", "x1")))
X[, 1] <- 1
y <- rpois(n, lambda = exp(X %*% c(0.5, 0.2)))

# Fit the Poisson mixture model
fit <- glmMixBayes(
  X = X, y = y, family = "poisson",
  control = list(iterations = 100, burnin.iterations = 50, seed = 602)
)
#> 
#> SAMPLING FOR MODEL 'glmMixBayes_poisson' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.34 seconds.
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
#> Chain 1: Iteration:  1 / 100 [  1%]  (Warmup)
#> Chain 1: Iteration: 10 / 100 [ 10%]  (Warmup)
#> Chain 1: Iteration: 20 / 100 [ 20%]  (Warmup)
#> Chain 1: Iteration: 30 / 100 [ 30%]  (Warmup)
#> Chain 1: Iteration: 40 / 100 [ 40%]  (Warmup)
#> Chain 1: Iteration: 50 / 100 [ 50%]  (Warmup)
#> Chain 1: Iteration: 51 / 100 [ 51%]  (Sampling)
#> Chain 1: Iteration: 60 / 100 [ 60%]  (Sampling)
#> Chain 1: Iteration: 70 / 100 [ 70%]  (Sampling)
#> Chain 1: Iteration: 80 / 100 [ 80%]  (Sampling)
#> Chain 1: Iteration: 90 / 100 [ 90%]  (Sampling)
#> Chain 1: Iteration: 100 / 100 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.049 seconds (Warm-up)
#> Chain 1:                0.13 seconds (Sampling)
#> Chain 1:                0.179 seconds (Total)
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
#>     . ECR-ITERATIVE-1                0.048                Converged (2 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

# Generate and print the comprehensive summary
fit_summary <- summary(fit)
print(fit_summary)
#> Call:
#> glmMixBayes(X = X, y = y, family = "poisson", control = list(iterations = 100, 
#>     burnin.iterations = 50, seed = 602))
#>  
#> Family:
#> poisson
#>  
#> (For Correct Matches):
#> Outcome Model Coefficients:
#>           Estimates Std. Error    2.5 % 97.5 %
#> Intercept   0.53209    0.07776  0.38256  0.697
#> x1          0.06383    0.06993 -0.04654  0.196
#>  
#> (For Mismatches):
#> Outcome Model Coefficients:
#>           Estimates Std. Error    2.5 % 97.5 %
#> Intercept   -0.5258     4.9104 -14.2191  8.715
#> x1           0.9850     5.7313  -9.7367 11.695
#>  
# }
```
