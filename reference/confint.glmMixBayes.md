# Credible intervals for regression coefficients from a glmMixBayes fit

Computes posterior credible intervals for the regression coefficients in
a fitted `glmMixBayes` model. By default, intervals are returned for all
coefficients in the primary component of the mixture model. A subset of
coefficients can be selected using `parm`.

## Usage

``` r
# S3 method for class 'glmMixBayes'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `glmMixBayes` model object.

- parm:

  Optional. Names or numeric indices of coefficients for which credible
  intervals should be returned. If `NULL`, intervals are returned for
  all coefficients.

- level:

  Probability level for the credible intervals. Defaults to `0.95`.

- ...:

  Not used.

## Value

A matrix with one row per coefficient and two columns giving the lower
and upper credible interval bounds. Row names correspond to coefficient
names.

## Examples

``` r
# Simulate simple Gaussian regression data
set.seed(604)
n <- 100

x1 <- rnorm(n)
X <- cbind(Intercept = 1, x1 = x1)

# True model: y = 1 + 2*x1 + error
y <- 1 + 2 * x1 + rnorm(n, sd = 1.5)

# Fit the model
fit <- glmMixBayes(
  X = X, y = y, family = "gaussian",
  control = list(iterations = 100, burnin.iterations = 50, seed = 604)
)
#> 
#> SAMPLING FOR MODEL 'glmMixBayes_gaussian' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 7.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.73 seconds.
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
#> Chain 1:  Elapsed Time: 0.075 seconds (Warm-up)
#> Chain 1:                0.079 seconds (Sampling)
#> Chain 1:                0.154 seconds (Total)
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
#>     . ECR-ITERATIVE-1                0.048                Converged (2 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

# Calculate 95% credible intervals for all coefficients
confint(fit, level = 0.95)
#>                2.5%    97.5%
#> Intercept 0.7260174 1.300971
#> x1        1.6197921 2.128374

# Extract the credible interval for the x1 coefficient
confint(fit, parm = "x1", level = 0.90)
#>         5%      95%
#> x1 1.67836 2.109774
```
