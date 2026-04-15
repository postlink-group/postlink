# Credible intervals for parameters from a survMixBayes fit

Computes posterior credible intervals for the regression coefficients,
mixing weight, and family-specific distribution parameters from a fitted
`survMixBayes` model. The returned intervals are organized by mixture
component and parameter type. Component 1 corresponds to the
correct-match component and component 2 to the incorrect-match
component.

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
component 2 is the incorrect-match component. Elements such as `theta`,
`shape1`, `shape2`, `scale1`, and `scale2` are numeric vectors of length
2 containing the lower and upper credible interval bounds for the
corresponding scalar parameters.

## Examples

``` r
# \donttest{
set.seed(301)
n <- 150
trt <- rbinom(n, 1, 0.5)

# Simulate Weibull AFT data
true_time <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))
cens_time <- rexp(n, rate = 0.1)
true_obs_time <- pmin(true_time, cens_time)
true_status <- as.integer(true_time <= cens_time)

# Induce linkage mismatch errors in approximately 20% of records
is_mismatch <- rbinom(n, 1, 0.2)
obs_time <- true_obs_time
obs_status <- true_status
mismatch_idx <- which(is_mismatch == 1)

shuffled <- sample(mismatch_idx)
obs_time[mismatch_idx] <- obs_time[shuffled]
obs_status[mismatch_idx] <- obs_status[shuffled]

linked_df <- data.frame(time = obs_time, status = obs_status, trt = trt)

adj <- adjMixBayes(linked.data = linked_df)

fit <- plsurvreg(
  survival::Surv(time, status) ~ trt,
  dist = "weibull",
  adjustment = adj,
  control = list(iterations = 200, burnin.iterations = 100, seed = 123)
)
#> 
#> SAMPLING FOR MODEL 'survMixBayes_weibull' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000125 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.25 seconds.
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
#> Chain 1:  Elapsed Time: 0.632 seconds (Warm-up)
#> Chain 1:                0.498 seconds (Sampling)
#> Chain 1:                1.13 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.18, indicating chains have not mixed.
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
#>     . ECR-ITERATIVE-1                0.078                Converged (2 iterations)       . 
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
#>            2.5%    97.5%
#> [1,] -0.7478373 1.989130
#> [2,]  0.2131905 1.253355
#> 
#> $coef2
#>           2.5%    97.5%
#> [1,] -3.145031 2.640371
#> [2,] -3.333294 2.923438
#> 
#> $theta
#>      2.5%     97.5% 
#> 0.5941939 0.9950833 
#> 
#> $shape1
#>      2.5%     97.5% 
#> 0.9991279 1.5116667 
#> 
#> $shape2
#>      2.5%     97.5% 
#> 0.6610367 4.9910123 
#> 
#> $scale1
#>      2.5%     97.5% 
#> 0.3702588 5.0051712 
#> 
#> $scale2
#>      2.5%     97.5% 
#> 0.5516955 5.0570164 
#> 

# Extract credible intervals specifically for the mixing weight
confint(fit, parm = "theta", level = 0.90)
#> $theta
#>        5%       95% 
#> 0.6192279 0.9942107 
#> 
# }
```
