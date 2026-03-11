# Credible intervals for glmMixBayes coefficients

Credible intervals for glmMixBayes coefficients

## Usage

``` r
# S3 method for class 'glmMixBayes'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `glmMixBayes` model object.

- parm:

  Optional. Parameter names or indices for selecting a subset of
  coefficients. If `NULL`, all coefficients are returned.

- level:

  Probability level for the intervals (default 0.95).

- ...:

  Not used.

## Value

A matrix with two columns (lower and upper bounds) and one row per
coefficient.

## Examples

``` r
# \donttest{
# Simulate data
set.seed(604)
n <- 100
X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("Intercept", "x1")))
X[, 1] <- 1
y <- rnorm(n, mean = X %*% c(1, 2), sd = 1.5)

# Fit the model
fit <- glmMixBayes(
  X = X, y = y, family = "gaussian",
  control = list(iterations = 100, burnin.iterations = 50, seed = 604)
)
#> 
#> SAMPLING FOR MODEL 'glmMixBayes_gaussian' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 7.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.75 seconds.
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
#> Chain 1:  Elapsed Time: 0.081 seconds (Warm-up)
#> Chain 1:                0.11 seconds (Sampling)
#> Chain 1:                0.191 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 2.19, indicating chains have not mixed.
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
#>     . ECR-ITERATIVE-1                0.05                 Converged (2 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

# Calculate 95% credible intervals for all outcome coefficients
confint(fit, level = 0.95)
#>                2.5%    97.5%
#> Intercept 0.7557846 1.346217
#> x1        1.5826242 2.093995

# Extract credible interval specifically for the 'x1' parameter
confint(fit, parm = "x1", level = 0.90)
#>          5%      95%
#> x1 1.613289 2.077088
# }
```
