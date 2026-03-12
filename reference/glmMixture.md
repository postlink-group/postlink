# GLM with Mixture-Based Linkage Error Adjustment

Fits a generalized linear model (GLM) accounting for mismatch errors
using a mixture model framework in the secondary analysis setting. The
variance-covariance matrix is estimated using the sandwich formula.

## Usage

``` r
glmMixture(
  x,
  y,
  family,
  z = cbind(rep(1, nrow(x))),
  m.rate = NULL,
  safe.matches = NULL,
  control = list(),
  ...
)
```

## Arguments

- x:

  Design matrix for the primary outcome model (numeric matrix or data
  frame).

- y:

  Response vector for the primary outcome model.

- family:

  A family object (e.g., `gaussian`, `binomial`) specifying the error
  distribution and link function. Can be a character string or a
  function.

- z:

  Design matrix for the mismatch indicator model (mismatch covariates).
  If NULL, an intercept-only model is assumed.

- m.rate:

  The assumed overall mismatch rate (a proportion between 0 and 1). If
  provided, it imposes a constraint on the mismatch model intercept.

- safe.matches:

  Logical vector; `TRUE` indicates a "safe match" (treated as definitely
  correct), `FALSE` indicates a potential mismatch.

- control:

  An optional list of control parameters. Arguments passed via `...`
  will override values in this list.

  - `max.iter`: Maximum EM iterations (default: 1000).

  - `cmax.iter`: Maximum iterations for the subroutine in the
    constrained logistic regression function (default: 1000).

  - `tol`: Convergence tolerance (default: 1e-4).

  - `init.beta`: Initial parameter estimates for the outcome model.

  - `init.gamma`: Initial parameter estimates for the mismatch indicator
    model.

  - `fy`: Estimated marginal density of the response. If NULL, estimated
    using Kernel Density Estimation or parametric assumption.

- ...:

  Additional arguments passed to `control`.

## Value

A list of results:

- coefficients:

  A named vector of coefficients for the outcome model.

- m.coefficients:

  A named vector of coefficients for the mismatch indicator model
  (gamma).

- match.prob:

  The posterior correct match probabilities (weights) for each
  observation.

- residuals:

  The working residuals, defined as `y - fitted.values`.

- fitted.values:

  The fitted mean values of the outcome model, obtained by transforming
  the linear predictors by the inverse of the link function.

- linear.predictors:

  The linear fit on the link scale.

- deviance:

  The deviance of the weighted outcome model at convergence.

- null.deviance:

  The deviance of the weighted null outcome model.

- var:

  The estimated variance-covariance matrix of the parameters (sandwich
  estimator).

- dispersion:

  The estimated dispersion parameter (e.g., variance for Gaussian,
  1/shape for Gamma).

- objective:

  A vector tracking the negative log pseudo-likelihood at each iteration
  of the EM algorithm.

- converged:

  Logical indicating if the EM algorithm converged within `max.iter`.

- rank:

  The numeric rank of the fitted linear model.

- df.residual:

  The residual degrees of freedom.

- df.null:

  The residual degrees of freedom for the null model.

- family:

  The `family` object used.

- `call`:

  The matched call.

## References

Slawski, M.\*, West, B. T., Bukke, P., Wang, Z., Diao, G., & Ben-David,
E. (2025). A general framework for regression with mismatched data based
on mixture modelling. *Journal of the Royal Statistical Society Series
A: Statistics in Society*, 188(3), 896-919.
[doi:10.1093/jrsssa/qnae083](https://doi.org/10.1093/jrsssa/qnae083)

Slawski, M.\*, Diao, G., Ben-David, E. (2021). A pseudo-likelihood
approach to linear regression with partially shuffled data. *Journal of
Computational and Graphical Statistics*. 30(4), 991-1003.
[doi:10.1080/10618600.2020.1870482](https://doi.org/10.1080/10618600.2020.1870482)

## Examples

``` r
data(lifem)

x <- cbind(1, poly(lifem$unit_yob, 3, raw = TRUE))
y <- lifem$age_at_death
z <- cbind(1, lifem$commf, lifem$comml)

fit <- glmMixture(x, y, family = "gaussian",
                  z, m.rate = 0.05, safe.matches = lifem$hndlnk)
```
