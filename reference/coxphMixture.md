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
#> Call:
#> NULL
#> 
#> Outcome Model Coefficients:
#>  d_cov  x_cov        
#>  7.005 -4.377 -2.899 
#> 
#> Mismatch Model Coefficients:
#> [1] 0.9644
#> 
#> Likelihood ratio test (model=outcome) not available due to pseudo-likelihood.
#> n= 200 , number of events= 144 
summary(fit)
#> 
#> Call:
#> NULL
#> 
#> --- Outcome Model (Cox PH) ---
#>             coef  exp(coef)   se(coef)       z Pr(>|z|)    
#> d_cov    7.00547 1102.65118    0.66757  10.494  < 2e-16 ***
#> x_cov   -4.37692    0.01256    0.37763 -11.590  < 2e-16 ***
#>         -2.89911    0.05507    0.61015  -4.751 2.02e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> --- Hazard Ratios & Confidence Intervals ---
#>       exp(coef) exp(-coef) lower .95 upper .95
#> d_cov 1.103e+03  9.069e-04 2.980e+02 4.080e+03
#> x_cov 1.256e-02  7.959e+01 5.994e-03 2.634e-02
#>       5.507e-02  1.816e+01 1.666e-02 1.821e-01
#> 
#> --- Mismatch Indicator Model ---
#>      Estimate Std. Error z value Pr(>|z|)    
#> [1,]   0.9644     0.1566   6.158 7.37e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Average Estimated Correct Match Rate: 0.8789 
#> Events: 144  / Total: 200 
#> Iterations: 24 
#> 
```
