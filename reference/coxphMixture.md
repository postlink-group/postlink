# CoxPH with Mixture-Based Linkage Error Adjustment

Fits a Cox proportional hazards regression adjusting for mismatched data
using a mixture modeling framework in the secondary analysis setting.
The method relies on a two-component mixture model where true matches
follow the Cox model and mismatches follow the marginal distribution of
the survival outcome. Variance estimates are obtained via Louis' method.

## Usage

``` r
coxphMixture(
  x,
  y,
  cens,
  z,
  m.rate = NULL,
  safe.matches = NULL,
  control = list(),
  ...
)
```

## Arguments

- x:

  A matrix or data.frame of covariates (design matrix).

- y:

  A numeric vector of observed time-to-event outcomes.

- cens:

  A numeric vector indicating censoring status (1 = censored, 0 =
  event). Note: This is the reverse of the standard `Surv` object
  convention where 1 usually indicates an event.

- z:

  A matrix or data.frame of mismatch covariates (e.g., match scores,
  blocking variables). Used to model the probability of a mismatch.

- m.rate:

  An optional numeric value between 0 and 1 specifying the assumed
  overall mismatch rate upper bound. If provided, the mismatch indicator
  model is constrained such that the average estimated mismatch rate
  does not exceed this bound.

- safe.matches:

  A logical vector indicating records known to be correct matches
  (TRUE). These records are fixed as matches (probability 1) during
  estimation. Defaults to all FALSE.

- control:

  An optional list of control parameters. Parameters can also be passed
  directly via `...`.

  - `louis.k`: Number of Monte Carlo iterations for variance estimation
    (default: 1000).

  - `max.iter`: Maximum EM iterations (default: 1000).

  - `cmax.iter`: Maximum iterations for the constrained optimization
    subroutine (default: 1000).

  - `tol`: Convergence tolerance (default: 1e-4).

  - `init.beta`: Initial estimates for outcome model coefficients.

  - `init.gamma`: Initial estimates for mismatch model coefficients.

  - `fy`: Pre-calculated marginal density of the response. If NULL,
    estimated non-parametrically.

- ...:

  Additional arguments passed to `control`.

## Value

An list of results:

- coefficients:

  Estimated coefficients for the outcome model (beta).

- m.coefficients:

  Estimated coefficients for the mismatch model (gamma).

- var:

  Variance-covariance matrix of the estimates.

- linear.predictors:

  Linear predictors for the outcome model.

- means:

  Column means of the covariate matrix `x`.

- n:

  Number of observations.

- nevent:

  Number of events.

- match.prob:

  Posterior probabilities that each observation is a correct match.

- objective:

  Value of the negative log pseudo-likelihood at each iteration.

- converged:

  Logical indicating if the algorithm converged.

- Lambdahat0:

  Estimated baseline cumulative hazard.

- gLambdahat0:

  the baseline cumulative hazard for the marginal density of the
  response variable (using Nelson-Aalen estimator)

## References

Bukke, P., Ben-David, E., Diao, G., Slawski, M., & West, B. T. (2025).
Cox Proportional Hazards Regression Using Linked Data: An Approach Based
on Mixture Modelling.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
set.seed(123)
n <- 200
# Generate covariates
x_cov <- seq(-3, 3, length = n)
d_cov <- rep(0:1, each = n/2)
X <- cbind(d_cov, x_cov, x_cov * d_cov)

# True parameters
b <- c(-1.5, 1, 0.5)
sigma <- 0.25
mu <- X %*% b
y <- exp(drop(mu)) * rweibull(n, shape = 1/sigma)

# Censoring
cens <- (y >= 1.5)
y[cens] <- 1.5

# Generate mismatch errors
ps <- rbeta(n, 4.5, 0.5)
logit_ps <- log(ps / (1 - ps))
mp <- cbind(1, logit_ps)
gamma_true <- c(-0.5, 1)
m <- 1 - rbinom(n, prob = plogis(mp %*% gamma_true), size = 1)
yperm <- y
shuffled_ix <- sample(which(m == 1))
yperm[shuffled_ix] <- yperm[sample(shuffled_ix)]

# Fit model
fit <- coxphMixture(x = X, y = yperm, cens = as.numeric(cens),
                    z = matrix(logit_ps, ncol = 1),
                    control = list(max.iter = 50))

print(fit)
summary(fit)
} # }
```
