# Print a brief summary of a glmMixBayes model

Print a brief summary of a glmMixBayes model

## Usage

``` r
# S3 method for class 'glmMixBayes'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `glmMixBayes`.

- digits:

  Minimum number of significant digits to show.

- ...:

  Further arguments (unused).

## Value

The input `x`, invisibly.

## Examples

``` r
# \donttest{
# 1. Simulate fast toy data (Gaussian mixture)
set.seed(601)
n <- 100
X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("Intercept", "x1")))
X[, 1] <- 1
y <- rnorm(n, mean = X %*% c(1, 2), sd = 1.5)

# 2. Fit the model with artificially low iterations for speed
fit <- glmMixBayes(
  X = X, y = y, family = "gaussian",
  control = list(iterations = 100, burnin.iterations = 50, seed = 601)
)
#> 
#> SAMPLING FOR MODEL 'glmMixBayes_gaussian' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.9e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.39 seconds.
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
#> Chain 1:  Elapsed Time: 0.058 seconds (Warm-up)
#> Chain 1:                0.102 seconds (Sampling)
#> Chain 1:                0.16 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 2.45, indicating chains have not mixed.
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
#>     . ECR-ITERATIVE-1                0.049                Converged (2 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

# 3. Print the model object
print(fit)
#> Call:
#> glmMixBayes(X = X, y = y, family = "gaussian", control = list(iterations = 100, 
#>     burnin.iterations = 50, seed = 601))
#> 
#> Coefficients:
#> Intercept         x1  
#>     1.038      2.218  
#> 
# }
```
