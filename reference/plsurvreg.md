# Fit Parametric Survival Models with Linkage Error Adjustment

`plsurvreg` fits parametric survival models (Accelerated Failure Time
models) to linked data. It serves as a wrapper for the internal engines,
compatible with
[`survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
specifications.

## Usage

``` r
plsurvreg(
  formula,
  adjustment,
  subset,
  na.action,
  dist = "weibull",
  model = TRUE,
  x = FALSE,
  y = FALSE,
  control = list(),
  ...
)
```

## Arguments

- formula:

  A formula object, with the response on the left of a ~ operator, and
  the terms on the right. The response must be a survival object as
  returned by the [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)
  function.

- adjustment:

  An object inheriting from class `"adjustment"`, or a `list` containing
  the necessary parameter specifications.

- subset:

  An optional vector specifying a subset of observations.

- na.action:

  A function for handling NAs.

- dist:

  Character string specifying the survival distribution (currently it
  must be "weibull" or "gamma").

- model:

  Logical; if `TRUE`, the model frame is returned.

- x, y:

  Logical; if `TRUE`, the model matrix (`x`) and response (`y`) are
  returned.

- control:

  A list of control parameters.

- ...:

  Additional arguments passed to the internal engine.

## Value

An object representing the fitted model. The specific class and
structure of the returned object depend directly on the `adjustment`
method provided:

- If `adjustment` is of class `adjMixBayes`, returns an object of class
  [`survregMixBayes`](https://postlink-group.github.io/postlink/reference/survregMixBayes.md).

## See also

[`adjMixBayes`](https://postlink-group.github.io/postlink/reference/adjMixBayes.md),
[`survregMixBayes`](https://postlink-group.github.io/postlink/reference/survregMixBayes.md)

## Examples

``` r
# \donttest{
library(survival)
set.seed(202)
n <- 200

# Simulate Weibull AFT data
trt <- rbinom(n, 1, 0.5)
true_time <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))
cens_time <- rexp(n, 0.05)
true_obs_time <- pmin(true_time, cens_time)
true_status <- as.numeric(true_time <= cens_time)

# Induce linkage mismatch errors by...
is_mismatch <- rbinom(n, 1, 0.2) # ~20% overall mismatch rate
obs_time <- true_obs_time
obs_status <- true_status
mismatch_idx <- which(is_mismatch == 1)

# Shuffle time and status together for mismatched records
shuffled <- sample(mismatch_idx)
obs_time[mismatch_idx] <- obs_time[shuffled]
obs_status[mismatch_idx] <- obs_status[shuffled]

linked_df <- data.frame(time = obs_time, status = obs_status, trt)

# Specify the Bayesian Mixture Adjustment
adj <- adjMixBayes(linked.data = linked_df)

# Fit the Adjusted Parametric Survival Model
fit <- plsurvreg(
  Surv(time, status) ~ trt,
  dist = "weibull",
  adjustment = adj,
  control = list(iterations = 2000, burnin.iterations = 500)
)
#> 
#> SAMPLING FOR MODEL 'survMixBayes_weibull' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000148 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.48 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  501 / 2000 [ 25%]  (Sampling)
#> Chain 1: Iteration:  700 / 2000 [ 35%]  (Sampling)
#> Chain 1: Iteration:  900 / 2000 [ 45%]  (Sampling)
#> Chain 1: Iteration: 1100 / 2000 [ 55%]  (Sampling)
#> Chain 1: Iteration: 1300 / 2000 [ 65%]  (Sampling)
#> Chain 1: Iteration: 1500 / 2000 [ 75%]  (Sampling)
#> Chain 1: Iteration: 1700 / 2000 [ 85%]  (Sampling)
#> Chain 1: Iteration: 1900 / 2000 [ 95%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 4.435 seconds (Warm-up)
#> Chain 1:                10.072 seconds (Sampling)
#> Chain 1:                14.507 seconds (Total)
#> Chain 1: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> 
#>     ......................................................................................
#>     . Method                         Time (sec)           Status                         . 
#>     ......................................................................................
#>     . ECR-ITERATIVE-1                1.363                Converged (3 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 1.4 seconds. 
# }
```
