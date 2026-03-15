# Credible intervals for parameters from a survMixBayes fit

Computes posterior credible intervals for the regression coefficients,
mixing weight, and family-specific distribution parameters from a fitted
`survMixBayes` model. The returned intervals are organized by mixture
component and parameter type. Component 1 corresponds to the
correct-match component and component 2 to the non-match component.

## Usage

``` r
# S3 method for class 'survMixBayes'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `survMixBayes`.

- parm:

  Optional character vector selecting which interval blocks to return.
  For example, `"theta"` returns only the credible interval for the
  mixing weight. If `NULL`, credible intervals are returned for all
  available parameter blocks.

- level:

  Probability level for the credible intervals. Defaults to `0.95`.

- ...:

  Further arguments (unused).

## Value

A named list of credible intervals. Elements `coef1` and `coef2` are
matrices with one row per regression coefficient and two columns giving
the lower and upper interval bounds for components 1 and 2,
respectively, where component 1 is the correct-match component and
component 2 is the non-match component. Elements such as `theta`,
`shape1`, `shape2`, `scale1`, and `scale2` are numeric vectors of length
2 containing the lower and upper credible interval bounds for the
corresponding scalar parameters.

## Examples

``` r
# \donttest{
# Simulate simple Weibull survival data
set.seed(403)
n <- 100

x1 <- rnorm(n)
x2 <- rnorm(n)
X <- cbind(x1 = x1, x2 = x2)

# Observed survival response: time and event indicator
y <- cbind(
  time = rweibull(n, shape = 1.2, scale = exp(x1)),
  event = rbinom(n, 1, 0.8)
)

# Fit the model
fit <- survregMixBayes(
  X = X, y = y, dist = "weibull",
  control = list(iterations = 100, burnin.iterations = 50, seed = 403)
)
#> 
#> SAMPLING FOR MODEL 'survMixBayes_weibull' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 7.9e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.79 seconds.
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
#> Chain 1:  Elapsed Time: 0.125 seconds (Warm-up)
#> Chain 1:                0.148 seconds (Sampling)
#> Chain 1:                0.273 seconds (Total)
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
#>     . ECR-ITERATIVE-1                0.047                Converged (2 iterations)       . 
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
