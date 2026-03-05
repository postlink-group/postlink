# Bayesian Two-Component Mixture Survival Regression (with label-switching adjustment)

Fit a two-component Bayesian parametric survival regression model in
Stan with component distributions `"gamma"` or `"weibull"` (both
components share the same family). Right-censored data are supported.

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

  A numeric design matrix (N x K), typically from
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  upstream. Missing values are not allowed.

- y:

  A survival response. Either a two-column numeric matrix with columns
  `time` and `event` (event indicator 1=event, 0=censored), or a list
  with elements `time` and `event`. Missing values are not allowed.

- dist:

  One of `"gamma"` or `"weibull"`. Controls the component-specific
  likelihood.

- priors:

  A named `list` (or `NULL`) of prior specifications. Because the Stan
  models are pre-compiled, these strings are parsed into numeric
  hyperparameters and passed to the model's data block. Any missing
  entries are automatically filled with symmetric defaults via
  `fill_defaults(priors, p_family = dist, model_type = "survival")`.

- control:

  A named `list` of tuning parameters with defaults:

  - `iterations` (default `1e4`) total iterations per chain;

  - `burnin.iterations` (default `1e3`) warm-up iterations;

  - `seed` (default random integer);

  - `cores` (default `getOption("mc.cores", 1L)`).

  Values in `...` override `control`.

- ...:

  Optional overrides for elements in `control`, e.g.
  `iterations = 4000`, `burnin.iterations = 1000`, `seed = 123`,
  `cores = 2`.

## Value

An object of class `"survMixBayes"` containing (at least):

- `m_samples`:

  Aligned \\z\\ label matrix (S x N).

- `estimates$coefficients`:

  Component 1 coefficient draws (S x K).

- `estimates$m.coefficients`:

  Component 2 coefficient draws (S x K).

- `estimates$theta`:

  Mixing weight draws for component 1 (length S).

- `estimates$shape`:

  Component-specific shape draws (family-specific).

- `estimates$m.shape`:

  Component-specific shape draws for component 2 (family-specific).

- `estimates$scale`:

  Component-specific scale draws (Weibull only).

- `estimates$m.scale`:

  Component-specific scale draws for component 2 (Weibull only).

- `family`:

  The survival family string.

- `call`:

  The matched call.

## Details

This implementation function assumes upstream wrappers have already
handled the formula/data interface and produced a design matrix `X` and
survival response `y`.

The function sets up the data, samples from the pre-compiled Stan
models, and then applies a label-switching correction (global majority
swap + ECR-ITERATIVE-1) to align component labels across MCMC draws.

## Label switching

We first perform an optional global swap \\(1 \leftrightarrow 2)\\ if
label 2 is more frequent overall, then align per-draw labels using
`ECR-ITERATIVE-1` permutations. Component-specific parameters are
permuted accordingly (e.g., `beta1`/`beta2`, `phi`/`shape`/`scale`, and
`theta`).

## Examples

``` r
# \donttest{
# 1. Simulate survival data from a two-component Weibull mixture
# Component 2 represents correct links (strong signal),
# and Component 1 represents mismatched links (noise).
set.seed(301)
n <- 150
X <- matrix(rnorm(n * 2), ncol = 2)
colnames(X) <- c("x1", "x2")

# Latent match status: 80% correct links (Z=2), 20% mismatches (Z=1)
Z_true <- rbinom(n, 1, 0.8) + 1

# Generate survival times based on latent status
time1 <- rweibull(n, shape = 1.2, scale = exp(0.1 * X[,1])) # Noise
time2 <- rweibull(n, shape = 1.5, scale = exp(0.5 * X[,1] - 0.5 * X[,2])) # Signal
obs_time <- ifelse(Z_true == 2, time2, time1)

# Apply right-censoring
cens_time <- rexp(n, rate = 0.1)
event <- as.integer(obs_time <= cens_time)
obs_time <- pmin(obs_time, cens_time)

y <- cbind(time = obs_time, event = event)

# 2. Fit the Bayesian Two-Component Mixture Survival Model
# Note: Iterations are set artificially low for run time.
fit <- survregMixBayes(
  X = X,
  y = y,
  dist = "weibull",
  control = list(iterations = 200, burnin.iterations = 100, seed = 123)
)
#> 
#> SAMPLING FOR MODEL 'survMixBayes_weibull' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 7.6e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.76 seconds.
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
#> Chain 1:  Elapsed Time: 0.33 seconds (Warm-up)
#> Chain 1:                0.323 seconds (Sampling)
#> Chain 1:                0.653 seconds (Total)
#> Chain 1: 
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
#>     . ECR-ITERATIVE-1                0.069                Converged (2 iterations)       . 
#>     ......................................................................................
#> 
#>     Relabelling all methods according to method ECR-ITERATIVE-1 ... done!
#>     Retrieve the 1 permutation arrays by typing:
#>         [...]$permutations$"ECR-ITERATIVE-1"
#>     Retrieve the 1 best clusterings: [...]$clusters
#>     Retrieve the 1 CPU times: [...]$timings
#>     Retrieve the 1 X 1 similarity matrix: [...]$similarity
#>     Label switching finished. Total time: 0.1 seconds. 

# 3. Inspect the aligned posterior estimates
# (Label switching is handled automatically via ECR-ITERATIVE-1)
cat("Component 2 (Correct Links) Coefficients:\n")
#> Component 2 (Correct Links) Coefficients:
print(colMeans(fit$estimates$m.coefficients))
#> [1] -0.5524656  0.6254333

cat("Component 1 (Mismatches) Coefficients:\n")
#> Component 1 (Mismatches) Coefficients:
print(colMeans(fit$estimates$coefficients))
#> [1]  0.2361357 -0.3872336

cat("Estimated mixing weight (Mismatch proportion):\n")
#> Estimated mixing weight (Mismatch proportion):
print(mean(fit$estimates$theta))
#> [1] 0.96202
# }
```
