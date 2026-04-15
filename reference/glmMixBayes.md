# Bayesian two-component mixture generalized linear model

Fits a Bayesian two-component mixture generalized linear model (GLM)
using Stan. Each observation is assumed to arise from one of two latent
components with component-specific regression coefficients.

## Usage

``` r
glmMixBayes(
  X,
  y,
  family = "gaussian",
  priors = NULL,
  control = list(iterations = 10000, burnin.iterations = 1000, seed =
    sample.int(.Machine$integer.max, 1), cores = getOption("mc.cores", 1L)),
  ...
)
```

## Arguments

- X:

  A numeric design matrix (N x K), typically from
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  upstream. Missing values are not allowed.

- y:

  A response vector of length N. For `"gaussian"` and `"gamma"`, `y`
  should be numeric; for `"poisson"` and `"binomial"`, `y` should be
  integer-valued (or coercible without loss).

- family:

  One of `"gaussian"`, `"poisson"`, `"binomial"`, or `"gamma"`. Controls
  the component-specific likelihood.

- priors:

  A named `list` (or `NULL`) of prior specifications. Because the Stan
  models are pre-compiled, these strings are parsed into numeric
  hyperparameters and passed to the model's data block. Any missing
  entries are automatically filled with symmetric defaults via
  `fill_defaults(priors, p_family = family, model_type = "glm")`.

- control:

  A named `list` of MCMC tuning parameters. Supported elements include:

  - `iterations`: total number of MCMC iterations per chain (default
    `1e4`);

  - `burnin.iterations`: number of warm-up iterations (default `1e3`);

  - `seed`: random seed used for reproducibility;

  - `cores`: number of CPU cores used for sampling (default
    `getOption("mc.cores", 1L)`).

  Values supplied through `...` override entries in `control`.

- ...:

  Optional arguments that override elements of `control`. For example,
  `iterations = 4000`, `burnin.iterations = 1000`, `seed = 123`, or
  `cores = 2`.

## Value

An object of class `"glmMixBayes"` containing (at least):

- `m_samples`:

  Posterior draws of aligned latent component labels (matrix of size
  draws × N), where component 1 corresponds to the correct-match
  component and component 2 to the incorrect-match component.

- `estimates$coefficients`:

  Posterior draws of regression coefficients for the correct-match
  component (component 1; draws × K).

- `estimates$m.coefficients`:

  Posterior draws of regression coefficients for the incorrect-match
  component (component 2; draws × K).

- `estimates$dispersion`:

  Posterior draws of the dispersion parameter for the correct-match
  component (component 1; family-specific).

- `estimates$m.dispersion`:

  Posterior draws of the dispersion parameter for the incorrect-match
  component (component 2; family-specific).

- `family`:

  The GLM family used in the model.

- `call`:

  The matched function call.

## Details

The function supports Gaussian, Poisson, Binomial, and Gamma outcome
families and returns posterior samples of the component-specific
regression parameters and mixture weight.

## Label switching

Mixture models are invariant to permutations of component labels, which
can lead to label switching in MCMC output. To ensure interpretable
posterior summaries, this function applies a post-processing step that
aligns component labels across posterior draws.

First, an optional global swap of labels (1 and 2) is performed if
component 2 is more frequent overall. Then, labels are aligned across
draws using the `ECR-ITERATIVE-1` relabeling algorithm.

## References

Gutman, R., Sammartino, C., Green, T., & Montague, B. (2016). Error
adjustments for file linking methods using encrypted unique client
identifier (eUCI) with application to recently released prisoners who
are HIV+. *Statistics in Medicine*, 35(1), 115–129.
[doi:10.1002/sim.6586](https://doi.org/10.1002/sim.6586)

Stephens, M. (2000). Dealing with label switching in mixture models.
*Journal of the Royal Statistical Society: Series B (Statistical
Methodology)*, 62(4), 795–809.
[doi:10.1111/1467-9868.00265](https://doi.org/10.1111/1467-9868.00265)

Papastamoulis, P. (2016). *label.switching*: An R package for dealing
with the label switching problem in MCMC outputs. *Journal of
Statistical Software*, 69(1), 1–24.
[doi:10.18637/jss.v069.c01](https://doi.org/10.18637/jss.v069.c01)

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

fit <- glmMixBayes(
  X = x,
  y = y,
  family = "gaussian",
  control = list(
    iterations = 200,
    burnin.iterations = 100,
    seed = 123
  )
)
#> 
#> SAMPLING FOR MODEL 'glmMixBayes_gaussian' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.52 seconds.
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
#> Chain 1:  Elapsed Time: 0.727 seconds (Warm-up)
#> Chain 1:                0.789 seconds (Sampling)
#> Chain 1:                1.516 seconds (Total)
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
#>     . ECR-ITERATIVE-1                0.07                 Converged (2 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 
# }
```
