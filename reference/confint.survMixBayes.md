# Confidence intervals for survMixBayes parameters

Returns posterior credible intervals for regression coefficients and key
parameters.

## Usage

``` r
# S3 method for class 'survMixBayes'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `survMixBayes`.

- parm:

  Optional. Character vector selecting which elements of the returned
  list to keep. If `NULL`, all intervals are returned.

- level:

  Confidence level.

- ...:

  Additional arguments (unused).

## Value

A named list of credible intervals. Elements `coef1` and `coef2` are
`p x 2` matrices (lower/upper) for component-specific regression
coefficients when available. Elements `theta`, `shape1`, `shape2`,
`scale1`, `scale2` are numeric vectors of length 2 giving lower/upper
credible intervals for scalar parameters.

## Examples

``` r
# \donttest{
# Simulate data
set.seed(403)
n <- 100
X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
y <- cbind(time = rweibull(n, shape = 1.2, scale = exp(X[, 1])),
           event = rbinom(n, 1, 0.8))

# Fit the model
fit <- survregMixBayes(
  X = X, y = y, dist = "weibull",
  control = list(iterations = 100, burnin.iterations = 50, seed = 403)
)
#> 
#> SAMPLING FOR MODEL 'survMixBayes_weibull' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 7.6e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.76 seconds.
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
#> Chain 1:  Elapsed Time: 0.124 seconds (Warm-up)
#> Chain 1:                0.147 seconds (Sampling)
#> Chain 1:                0.271 seconds (Total)
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

# Calculate 95% credible intervals for all parameters
confint(fit, level = 0.95)
#> $coef1
#>            2.5%     97.5%
#> [1,]  0.7937272 1.2774170
#> [2,] -0.2868521 0.3001631
#> 
#> $coef2
#>           2.5%    97.5%
#> [1,] -3.620008 4.379870
#> [2,] -2.081506 3.310001
#> 
#> $theta
#>      2.5%     97.5% 
#> 0.7056352 0.9997981 
#> 
#> $shape1
#>      2.5%     97.5% 
#> 0.9006904 1.2965417 
#> 
#> $shape2
#>      2.5%     97.5% 
#> 0.3646252 5.0810257 
#> 
#> $scale1
#>      2.5%     97.5% 
#> 0.8174115 1.4282688 
#> 
#> $scale2
#>       2.5%      97.5% 
#> 0.08103639 4.10832947 
#> 

# Extract credible intervals specifically for the mixing weight
confint(fit, parm = "theta", level = 0.90)
#> $theta
#>        5%       95% 
#> 0.7253824 0.9997456 
#> 
# }
```
