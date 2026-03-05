# Multiple Imputation Pooling from Component Allocation Posterior Samples

Pool regression estimates across many "completed datasets" that arise
when each posterior draw of component allocation indicators \\z\\
selects a subset of records to treat as component 1. For each retained
draw, the function fits the requested model (LM/GLM) on the selected
subset, collects coefficients and their covariance, and combines them
using Rubin's rules (vector/matrix form).

## Usage

``` r
# S3 method for class 'glmMixBayes'
mi_with(
  object,
  data,
  formula,
  family = NULL,
  min_n = NULL,
  quietly = TRUE,
  ...
)
```

## Arguments

- object:

  A `glmMixBayes` model object containing posterior allocation samples.

- data:

  A data.frame with all candidate records in the same row order as used
  in the model.

- formula:

  Model formula for refitting on each draw (required).

- family:

  A [`stats::family()`](https://rdrr.io/r/stats/family.html) object;
  defaults to a canonical family derived from `object$family`.

- min_n:

  Minimum sample size required to fit the model for a given draw.
  Defaults to `p + 1`, where `p` is the number of columns in the model
  matrix.

- quietly:

  If `TRUE`, suppress errors from individual failed fits and skip them.

- ...:

  Additional arguments passed through (currently unused).

## Value

An object of class `c("mi_link_pool_glm", "mi_link_pool")`.

## Examples

``` r
# \donttest{
# 1. Simulate data linked with errors
set.seed(606)
n <- 100
linked_data <- data.frame(
  x1 = rnorm(n),
  y = rnorm(n)
)
X <- model.matrix(~ x1, data = linked_data)

# 2. Fit the GLM mixture model
fit <- glmMixBayes(
  X = X, y = linked_data$y, family = "gaussian",
  control = list(iterations = 150, burnin.iterations = 50, seed = 606)
)
#> 
#> SAMPLING FOR MODEL 'glmMixBayes_gaussian' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.4 seconds.
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
#> Chain 1: Iteration:   1 / 150 [  0%]  (Warmup)
#> Chain 1: Iteration:  15 / 150 [ 10%]  (Warmup)
#> Chain 1: Iteration:  30 / 150 [ 20%]  (Warmup)
#> Chain 1: Iteration:  45 / 150 [ 30%]  (Warmup)
#> Chain 1: Iteration:  51 / 150 [ 34%]  (Sampling)
#> Chain 1: Iteration:  65 / 150 [ 43%]  (Sampling)
#> Chain 1: Iteration:  80 / 150 [ 53%]  (Sampling)
#> Chain 1: Iteration:  95 / 150 [ 63%]  (Sampling)
#> Chain 1: Iteration: 110 / 150 [ 73%]  (Sampling)
#> Chain 1: Iteration: 125 / 150 [ 83%]  (Sampling)
#> Chain 1: Iteration: 140 / 150 [ 93%]  (Sampling)
#> Chain 1: Iteration: 150 / 150 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.065 seconds (Warm-up)
#> Chain 1:                0.211 seconds (Sampling)
#> Chain 1:                0.276 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.15, indicating chains have not mixed.
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
#>     . ECR-ITERATIVE-1                0.072                Converged (2 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

# 3. Multiple Imputation Pooling
# For each MCMC draw, the model extracts the latent correct-match indicators,
# fits the specified GLM on the subset of "correct" matches, and pools the
# results using Rubin's rules.
pooled_fit <- mi_with(
  object = fit,
  data = linked_data,
  formula = y ~ x1,
  family = gaussian()
)

# 4. View pooled results
print(pooled_fit)
#> Multiple Imputation (from allocation draws):
#>   Retained imputations (m): 100 
#> 
#>             Estimate Std.Error   CI.lwr  CI.upr       df
#> (Intercept)  0.16500   0.18639 -0.20212 0.53213 245.6075
#> x1           0.03584   0.15630 -0.27159 0.34328 340.6614
# }
```
