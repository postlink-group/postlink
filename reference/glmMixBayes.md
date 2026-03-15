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
if (FALSE) { # \dontrun{
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
} # }
```
