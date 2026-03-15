# Summary method for survMixBayes models

Computes posterior summaries for the regression coefficients, mixing
weight, and component-specific distribution parameters in a fitted
`survMixBayes` model. Throughout, component 1 is interpreted as the
correct-match component and component 2 as the non-match component.

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
correct-match component and component 2 to the non-match component.

## Examples

``` r
# \donttest{
# Simulate data
set.seed(402)
n <- 100
X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
y <- cbind(time = rweibull(n, shape = 1.2, scale = exp(X[, 1])),
           event = rbinom(n, 1, 0.8))

# Fit the model
fit <- survregMixBayes(
  X = X, y = y, dist = "weibull",
  control = list(iterations = 100, burnin.iterations = 50, seed = 402)
)
#> 
#> SAMPLING FOR MODEL 'survMixBayes_weibull' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 6.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.65 seconds.
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
#> Chain 1:  Elapsed Time: 0.059 seconds (Warm-up)
#> Chain 1:                0.062 seconds (Sampling)
#> Chain 1:                0.121 seconds (Total)
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
#>     . ECR-ITERATIVE-1                0.097                Converged (4 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

# Generate and print posterior summaries
fit_summary <- summary(fit, probs = c(0.025, 0.5, 0.975))
print(fit_summary)
#> Summary of Bayesian mixture survival regression
#> 
#> Call:
#> survregMixBayes(X = X, y = y, dist = "weibull", control = list(iterations = 100, 
#>     burnin.iterations = 50, seed = 402))
#> 
#> Family: weibull 
#> 
#> Coefficients (component 1 = correct-match) (quantiles):
#>         2.5%    50%  97.5%
#> [1,]  0.8265 1.0880 1.4977
#> [2,] -0.3615 0.0267 0.5531
#> 
#> Coefficients (component 2 = incorrect-match) (quantiles):
#>         2.5%    50%  97.5%
#> [1,]  0.0463 0.9834 1.5892
#> [2,] -0.5977 0.4045 1.0848
#> 
#> Theta (mix weight for component 1 = correct-match):
#>   2.5%    50%  97.5% 
#> 0.5240 0.7437 0.9064 
#> 
#> Shape (component 1 = correct-match) (quantiles):
#>        2.5%    50%  97.5%
#> [1,] 1.1509 1.4914 2.0428
#> 
#> Shape (component 2 = incorrect-match) (quantiles):
#>        2.5%    50%  97.5%
#> [1,] 0.9844 1.6636 5.1042
#> 
#> Scale (component 1 = correct-match) (quantiles):
#>       2.5%    50%  97.5%
#> [1,] 0.889 1.2273 1.6315
#> 
#> Scale (component 2 = incorrect-match) (quantiles):
#>        2.5%    50%  97.5%
#> [1,] 0.5543 1.1801 3.0561
# }
```
