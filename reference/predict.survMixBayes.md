# Predictions from a survMixBayes model

Computes posterior predictions for each latent component of a
`survMixBayes` model. By default, predictions are returned on the linear
predictor scale for both components.

## Usage

``` r
# S3 method for class 'survMixBayes'
predict(
  object,
  newdata = NULL,
  se.fit = FALSE,
  interval = c("none", "credible"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `survMixBayes` model object.

- newdata:

  A numeric matrix of new observations (\\n\_{new} \times K\\) with
  columns aligned to the design matrix used for fitting. If `NULL`, the
  fitted design matrix stored in `object$X` is used.

- se.fit:

  Logical; if `TRUE`, also return posterior SD of predictions.

- interval:

  Either `"none"` or `"credible"`, indicating whether to compute
  credible intervals.

- level:

  Probability level for the credible interval (default 0.95).

- ...:

  Not used.

## Value

A list with two components, `component1` and `component2`, corresponding
to the two latent mixture components. If `se.fit = FALSE` and
`interval = "none"`, each element is a numeric vector of posterior mean
linear predictors. Otherwise, each element is a matrix containing the
fitted values and, optionally, posterior SDs and credible interval
bounds.

## Details

Component 1 is interpreted as the correct-match component and component
2 as the incorrect-match component (after label-switching correction).

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
  control = list(
    iterations = 200,
    burnin.iterations = 100,
    seed = 123
  )
)
#> 
#> SAMPLING FOR MODEL 'survMixBayes_weibull' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000104 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.04 seconds.
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
#> Chain 1:  Elapsed Time: 0.642 seconds (Warm-up)
#> Chain 1:                0.505 seconds (Sampling)
#> Chain 1:                1.147 seconds (Total)
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

# Create a design matrix for new covariate values
newdata <- stats::model.matrix(~ trt, data = data.frame(trt = c(0, 1)))

# Predict posterior mean linear predictors for each latent component
preds <- predict(fit, newdata = newdata, se.fit = TRUE, interval = "credible")
print(preds$component1)
#>         fit    se.fit       2.5 %   97.5 %
#> 1 0.4407128 0.7432820 -0.74783735 1.989130
#> 2 1.1569081 0.7328594 -0.04520797 2.696848
print(preds$component2)
#>          fit   se.fit     2.5 %   97.5 %
#> 1 -0.0683974 1.550691 -3.145031 2.640371
#> 2 -0.2061190 2.156831 -4.552259 2.747006
# }
```
