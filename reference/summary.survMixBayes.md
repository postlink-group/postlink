# Summary method for survMixBayes models

Computes posterior summaries for the regression coefficients, mixing
weight, and component-specific distribution parameters in a fitted
`survMixBayes` model. Throughout, component 1 is interpreted as the
correct-match component and component 2 as the incorrect-match
component.

## Usage

``` r
# S3 method for class 'survMixBayes'
summary(object, probs = c(0.025, 0.5, 0.975), ...)
```

## Arguments

- object:

  An object of class `survMixBayes`.

- probs:

  Numeric vector of probabilities used to compute posterior quantiles
  for the model parameters. The default, `c(0.025, 0.5, 0.975)`, gives a
  posterior median and a 95\\ credible interval.

- ...:

  Further arguments (unused).

## Value

An object of class `summary.survMixBayes` containing posterior quantile
summaries for the regression coefficients in both mixture components,
the mixing weight, and any family-specific distribution parameters
included in the fitted model. Component 1 corresponds to the
correct-match component and component 2 to the incorrect-match
component.

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
#> Chain 1: Gradient evaluation took 9.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.95 seconds.
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
#> Chain 1:  Elapsed Time: 0.618 seconds (Warm-up)
#> Chain 1:                0.488 seconds (Sampling)
#> Chain 1:                1.106 seconds (Total)
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
#>     . ECR-ITERATIVE-1                0.104                Converged (2 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

fit_summary <- summary(fit, probs = c(0.025, 0.5, 0.975))
print(fit_summary)
#> Summary of Bayesian mixture survival regression
#> 
#> Call:
#> plsurvreg(formula = survival::Surv(time, status) ~ trt, adjustment = adj, 
#>     dist = "weibull", control = list(iterations = 200, burnin.iterations = 100, 
#>         seed = 123))
#> 
#> Family: weibull 
#> 
#> Coefficients (component 1 = correct-match) (quantiles):
#>         2.5%    50%  97.5%
#> [1,] -0.7478 0.3784 1.9891
#> [2,]  0.2132 0.7247 1.2534
#> 
#> Coefficients (component 2 = incorrect-match) (quantiles):
#>         2.5%     50%  97.5%
#> [1,] -3.1450 -0.0303 2.6404
#> [2,] -3.3333  0.1177 2.9234
#> 
#> Theta (mix weight for component 1 = correct-match):
#>   2.5%    50%  97.5% 
#> 0.5942 0.9379 0.9951 
#> 
#> Shape (component 1 = correct-match) (quantiles):
#>        2.5%    50%  97.5%
#> [1,] 0.9991 1.2636 1.5117
#> 
#> Shape (component 2 = incorrect-match) (quantiles):
#>       2.5%    50% 97.5%
#> [1,] 0.661 1.6967 4.991
#> 
#> Scale (component 1 = correct-match) (quantiles):
#>        2.5%    50%  97.5%
#> [1,] 0.3703 1.9828 5.0052
#> 
#> Scale (component 2 = incorrect-match) (quantiles):
#>        2.5%    50% 97.5%
#> [1,] 0.5517 1.9333 5.057
# }
```
