# Predictions from a glmMixBayes model

Predictions from a glmMixBayes model

## Usage

``` r
# S3 method for class 'glmMixBayes'
predict(
  object,
  newx,
  type = c("link", "response"),
  se.fit = FALSE,
  interval = c("none", "credible"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `glmMixBayes` model object.

- newx:

  A numeric matrix of new observations (n_new x K) with columns aligned
  to the design matrix `X` used for fitting.

- type:

  Either `"link"` or `"response"`, indicating the scale of predictions.

- se.fit:

  Logical; if `TRUE`, also return posterior SD of predictions.

- interval:

  Either `"none"` or `"credible"`, indicating whether to compute a
  credible interval.

- level:

  Probability level for the credible interval (default 0.95).

- ...:

  Not used.

## Value

If `se.fit = FALSE` and `interval = "none"`, a numeric vector of
predicted values. Otherwise, a matrix with columns for the fit,
(optional) `se.fit`, and (optional) credible interval bounds.

## Examples

``` r
# \donttest{
# Simulate data
set.seed(605)
n <- 100
X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("Intercept", "x1")))
X[, 1] <- 1
y <- rpois(n, lambda = exp(X %*% c(0.5, 0.8)))

# Fit the Poisson mixture model
fit <- glmMixBayes(
  X = X, y = y, family = "poisson",
  control = list(iterations = 100, burnin.iterations = 50, seed = 605)
)
#> 
#> SAMPLING FOR MODEL 'glmMixBayes_poisson' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.51 seconds.
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
#> Chain 1:  Elapsed Time: 0.04 seconds (Warm-up)
#> Chain 1:                0.204 seconds (Sampling)
#> Chain 1:                0.244 seconds (Total)
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
#>     . ECR-ITERATIVE-1                0.037                Converged (2 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0 seconds. 

# Create a new design matrix for prediction
X_new <- matrix(c(1, 1,   # Intercepts
                  0, 1.5), # x1 values
                ncol = 2, dimnames = list(NULL, c("Intercept", "x1")))

# Predict expected counts (type = "response") with 95% credible intervals
preds <- predict(fit, newx = X_new, type = "response", interval = "credible")
print(preds)
#>           fit    2.5 %   97.5 %
#> [1,] 1.540914 1.289699 1.817937
#> [2,] 4.864849 4.367653 5.438319
# }
```
