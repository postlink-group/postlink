# Pool posterior component classifications from a survMixBayes fit

Summarizes posterior component classifications from a fitted
`survMixBayes` model by computing, for each observation, the posterior
probability of belonging to component 1. In this package, component 1 is
interpreted as the correct-match component.

## Usage

``` r
# S3 method for class 'survMixBayes'
mi_with(object, ...)
```

## Arguments

- object:

  An object of class `survMixBayes`.

- ...:

  Further arguments (unused).

## Value

An object of class `c("mi_link_pool_survreg", "mi_link_pool")`
containing the fitted model call, survival distribution, posterior
probabilities of belonging to component 1 for each observation, where
component 1 is interpreted as the correct-match component, and the
stored posterior draws from the original fitted model.

## Examples

``` r
# \donttest{
# Simulate data
set.seed(406)
n <- 100
X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
y <- cbind(time = rweibull(n, shape = 1.2, scale = exp(X[, 1])),
           event = rbinom(n, 1, 0.8))

# Fit the model
fit <- survregMixBayes(
  X = X, y = y, dist = "weibull",
  control = list(iterations = 100, burnin.iterations = 50, seed = 406)
)
#> 
#> SAMPLING FOR MODEL 'survMixBayes_weibull' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 6.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.63 seconds.
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
#> Chain 1:  Elapsed Time: 0.062 seconds (Warm-up)
#> Chain 1:                0.094 seconds (Sampling)
#> Chain 1:                0.156 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.28, indicating chains have not mixed.
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

# Summarize posterior component classifications
pooled_obj <- mi_with(fit)

# View the posterior probabilities of belonging to component 1
summary(pooled_obj$p_component1)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   0.340   0.820   0.880   0.847   0.920   1.000 
# }
```
