# Bayesian Two-Component Mixture GLM (with label-switching adjustment)

Fit a two-component Bayesian mixture generalized linear model (GLM) in
Stan with families `"gaussian"`, `"poisson"`, `"binomial"`, or
`"gamma"`. This implementation function assumes upstream wrappers have
already handled the formula/data interface and produced a design matrix
`X` and response `y`.

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

An object of class `"glmMixBayes"` containing (at least):

- `m_samples`:

  Aligned \\z\\ label matrix (S x N).

- `estimates$coefficients`:

  Component 1 coefficient draws (S x K).

- `estimates$m.coefficients`:

  Component 2 coefficient draws (S x K).

- `estimates$dispersion`:

  Component 1 dispersion (family-specific).

- `estimates$m.dispersion`:

  Component 2 dispersion (family-specific).

- `family`:

  The GLM family string.

- `call`:

  The matched call.

## Details

The function builds Stan code, compiles, samples, and then applies a
label-switching correction (global majority swap + ECR-ITERATIVE-1) to
align component labels across MCMC draws.

## Label switching

We first perform an optional global swap \\(1 \leftrightarrow 2)\\ if
label 2 is more frequent overall, then align per-draw labels using
`ECR-ITERATIVE-1` permutations. Component-specific parameters are
permuted accordingly (e.g., `beta1`/`beta2`, and `sigma`/`phi`).

## Examples

``` r
# \donttest{
# 1. Simulate data from a two-component Gaussian mixture
# Component 2: Correct links (strong signal/relationship)
# Component 1: False links (noise/null relationship)
set.seed(501)
n <- 150
X <- matrix(rnorm(n * 2), ncol = 2)
colnames(X) <- c("Intercept", "x1")
X[, 1] <- 1 # Set intercept

# Latent match status: 70% correct links (Z=2), 30% mismatches (Z=1)
Z_true <- rbinom(n, 1, 0.7) + 1

# Generate responses based on latent status
y1 <- rnorm(n, mean = X %*% c(0, 0), sd = 2)     # Noise (mismatches)
y2 <- rnorm(n, mean = X %*% c(1, 1.5), sd = 0.5) # Signal (correct links)
y <- ifelse(Z_true == 2, y2, y1)

# 2. Fit the Bayesian Two-Component Mixture GLM
# Note: Iterations are set artificially low for check speed.
fit <- glmMixBayes(
  X = X,
  y = y,
  family = "gaussian",
  control = list(iterations = 200, burnin.iterations = 100, seed = 123)
)
#> 
#> SAMPLING FOR MODEL 'glmMixBayes_gaussian' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.57 seconds.
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
#> Chain 1:  Elapsed Time: 0.089 seconds (Warm-up)
#> Chain 1:                0.07 seconds (Sampling)
#> Chain 1:                0.159 seconds (Total)
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
#> Global label swap performed: label 2 dominates label 1.
#> 
#>     ......................................................................................
#>     . Method                         Time (sec)           Status                         . 
#>     ......................................................................................
#>     . ECR-ITERATIVE-1                0.087                Converged (2 iterations)       . 
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
#>  Intercept         x1 
#>  0.3403996 -1.2393847 

cat("Component 1 (Mismatches) Coefficients:\n")
#> Component 1 (Mismatches) Coefficients:
print(colMeans(fit$estimates$coefficients))
#> Intercept        x1 
#>  1.049480  1.566094 
# }
```
