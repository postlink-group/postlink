# Pooling regression fits across posterior draws of correct-match classifications

Use posterior draws of the latent match indicators from
[`glmMixBayes()`](https://postlink-group.github.io/postlink/reference/glmMixBayes.md)
to repeatedly identify which records are treated as correct matches,
refit the requested regression model on those records, and pool the
resulting estimates.

Each retained posterior draw defines one subset of records classified as
correct matches. The function fits the specified
[`lm()`](https://rdrr.io/r/stats/lm.html) or
[`glm()`](https://rdrr.io/r/stats/glm.html) model to that subset,
extracts the estimated coefficients and their covariance matrix, and
combines the results across draws using multiple-imputation pooling
rules.

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

  A `glmMixBayes` model object containing posterior draws of the latent
  match indicators.

- data:

  A data.frame with all candidate records in the same row order as used
  in the model.

- formula:

  Model formula for refitting on each draw (required).

- family:

  A [`stats::family()`](https://rdrr.io/r/stats/family.html) object for
  the refitted model. If not supplied, the function chooses a default
  family based on `object$family`.

- min_n:

  Minimum number of records required to fit the model for a given
  posterior draw. The default is `p + 1`, where `p` is the number of
  columns in the model matrix.

- quietly:

  If `TRUE`, draws that lead to fitting errors are skipped without
  printing the full error message.

- ...:

  Additional arguments passed through (currently unused).

## Value

An object of class `c("mi_link_pool_glm", "mi_link_pool")` containing
pooled coefficient estimates, standard errors, confidence intervals, and
related summary information.

## Examples

``` r
# \donttest{
data(lifem)

# lifem data preprocessing
# For computational efficiency in the example, we work with a subset of the lifem data.
lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
lifem_small <- rbind(
  head(subset(lifem, hndlnk == 1), 100),
  head(subset(lifem, hndlnk == 0), 20)
)

x <- cbind(1, poly(lifem_small$unit_yob, 3, raw = TRUE))
y <- lifem_small$age_at_death

adj <- adjMixBayes(
  linked.data = lifem_small,
  priors = list(theta = "beta(2, 2)")
)

fit <- plglm(
  age_at_death ~ poly(unit_yob, 3, raw = TRUE),
  family = "gaussian",
  adjustment = adj,
  control = list(
    iterations = 200,
    burnin.iterations = 100,
    seed = 123
  )
)
#> 
#> SAMPLING FOR MODEL 'glmMixBayes_gaussian' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.9e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.49 seconds.
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
#> Chain 1:  Elapsed Time: 0.599 seconds (Warm-up)
#> Chain 1:                0.769 seconds (Sampling)
#> Chain 1:                1.368 seconds (Total)
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

pooled_fit <- mi_with(
  object = fit,
  data = lifem_small,
  formula = age_at_death ~ poly(unit_yob, 3, raw = TRUE),
  family = gaussian()
)

print(pooled_fit)
#> Pooled regression results across posterior match classifications:
#>   Retained imputations (m): 100 
#> 
#>                                 Estimate Std.Error     CI.lwr    CI.upr
#> (Intercept)                     59.84352   4.14959   51.71046  67.97658
#> poly(unit_yob, 3, raw = TRUE)1 -21.86377  41.24687 -102.70915  58.98161
#> poly(unit_yob, 3, raw = TRUE)2 -18.42207 108.61197 -231.31841 194.47428
#> poly(unit_yob, 3, raw = TRUE)3  64.77724  75.61585  -83.44377 212.99825
#>                                        df
#> (Intercept)                    1127335.70
#> poly(unit_yob, 3, raw = TRUE)1   32513.19
#> poly(unit_yob, 3, raw = TRUE)2   12392.49
#> poly(unit_yob, 3, raw = TRUE)3   10769.54
# }
```
