# Bayesian two-component mixture survival regression model

Fits a Bayesian two-component parametric survival regression model using
Stan. Each observation is assumed to arise from one of two latent
components with component-specific survival regression parameters.

## Usage

``` r
survregMixBayes(
  X,
  y,
  dist = "weibull",
  priors = NULL,
  control = list(iterations = 10000, burnin.iterations = 1000, seed =
    sample.int(.Machine$integer.max, 1), cores = getOption("mc.cores", 1L)),
  ...
)
```

## Arguments

- X:

  A numeric design matrix (\\N \times K\\), typically created by
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html). Each
  row corresponds to one observation and each column to one covariate in
  the survival model. Missing values are not allowed.

- y:

  A survival response. This can be either a two-column numeric matrix
  with columns `time` and `event`, where `event = 1` indicates an
  observed event and `event = 0` indicates right censoring, or a list
  with elements `time` and `event`. Missing values are not allowed.

- dist:

  Character string specifying the parametric survival distribution used
  for both mixture components. Supported values are `"gamma"` and
  `"weibull"`.

- priors:

  A named `list` of prior specifications, or `NULL`. Since the Stan
  models are pre-compiled, prior specifications are converted into the
  corresponding numeric hyperparameters and passed to the model as data.
  Any missing entries are automatically filled in using symmetric
  default values via
  `fill_defaults(priors, p_family = dist, model_type = "survival")`.

- control:

  A named `list` of control parameters for posterior sampling. Defaults
  are:

  - `iterations` (default `1e4`): total number of iterations per chain;

  - `burnin.iterations` (default `1e3`): number of warm-up iterations;

  - `seed` (default: a random integer): random seed for reproducibility;

  - `cores` (default `getOption("mc.cores", 1L)`): number of CPU cores
    used.

  Values supplied through `...` override the corresponding entries in
  `control`.

- ...:

  Optional overrides for elements of `control`, such as
  `iterations = 4000`, `burnin.iterations = 1000`, `seed = 123`, or
  `cores = 2`.

## Value

An object of class `"survMixBayes"` containing (at least):

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

- `estimates$theta`:

  Posterior draws of the mixing weight for the correct-match component
  (component 1; vector of length draws).

- `estimates$shape`:

  Posterior draws of the shape parameter for the correct-match component
  (component 1; family-specific).

- `estimates$m.shape`:

  Posterior draws of the shape parameter for the incorrect-match
  component (component 2; family-specific).

- `estimates$scale`:

  Posterior draws of the scale parameter for the correct-match component
  (component 1; Weibull only).

- `estimates$m.scale`:

  Posterior draws of the scale parameter for the incorrect-match
  component (component 2; Weibull only).

- `family`:

  The survival distribution used in the model.

- `call`:

  The matched function call.

## Details

The function supports `"gamma"` and `"weibull"` component distributions,
with both components sharing the same family. Right-censored survival
outcomes are supported.

Posterior draws are returned for the component-specific regression
parameters and mixing weight. To improve interpretability of posterior
summaries, the function applies a post-processing step that aligns
component labels across posterior draws.

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
# Example: Bayesian mixture survival model fit to linked survival data
# with induced linkage mismatch errors

# 1. Simulate a linked survival dataset
set.seed(301)
n <- 150
X <- matrix(rnorm(n * 2), ncol = 2)
colnames(X) <- c("x1", "x2")

# Generate survival times from a Weibull AFT model
true_time <- rweibull(
  n,
  shape = 1.5,
  scale = exp(0.5 * X[, 1] - 0.5 * X[, 2])
)

# Apply right-censoring
cens_time <- rexp(n, rate = 0.1)
event <- as.integer(true_time <= cens_time)
obs_time <- pmin(true_time, cens_time)

# Induce linkage mismatch errors in approximately 20% of records
is_mismatch <- rbinom(n, 1, 0.2)
mismatch_idx <- which(is_mismatch == 1)

shuffled <- sample(mismatch_idx)
obs_time[mismatch_idx] <- obs_time[shuffled]
event[mismatch_idx] <- event[shuffled]

y <- cbind(time = obs_time, event = event)

# 2. Fit the Bayesian two-component mixture survival model
# Note: Iterations are set artificially low for run time
fit <- survregMixBayes(
  X = X,
  y = y,
  dist = "weibull",
  control = list(iterations = 200, burnin.iterations = 100, seed = 123)
)
#> 
#> SAMPLING FOR MODEL 'survMixBayes_weibull' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 9.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.94 seconds.
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
#> Chain 1:  Elapsed Time: 0.194 seconds (Warm-up)
#> Chain 1:                0.185 seconds (Sampling)
#> Chain 1:                0.379 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.19, indicating chains have not mixed.
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
#>     . ECR-ITERATIVE-1                0.134                Converged (3 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

# 3. Inspect posterior summaries
# (Label switching is handled automatically)
cat("Component 1 (Correct Links):\n")
#> Component 1 (Correct Links):
print(colMeans(fit$estimates$coefficients))
#> [1]  0.3856700 -0.3266702

cat("Component 2 (Incorrect Links):\n")
#> Component 2 (Incorrect Links):
print(colMeans(fit$estimates$m.coefficients))
#> [1]  0.1048798 -0.1195878

cat("Estimated probability of correct linkage:\n")
#> Estimated probability of correct linkage:
print(mean(fit$estimates$theta))
#> [1] 0.8279494
# }
```
