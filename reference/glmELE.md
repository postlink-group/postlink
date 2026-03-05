# GLM with ELE-Based Linkage Error Adjustment

Fits a generalized linear model (GLM) accounting for exchangeable
linkage errors (ELE) as defined by Chambers (2009). It solves the
unbiased estimating equations resulting from the modified mean function
induced by mismatch errors.

## Usage

``` r
glmELE(
  x,
  y,
  family = "gaussian",
  m.rate,
  audit.size = NULL,
  blocks,
  weight.matrix = "all",
  control = list(init.beta = NULL),
  ...
)
```

## Arguments

- x:

  A numeric matrix of predictors (design matrix).

- y:

  A numeric vector of responses.

- family:

  the type of regression model ("gaussian" - default, "poisson",
  "binomial", "gamma"). Standard link functions are used ("identity" for
  Gaussian, "log" for Poisson and Gamma, and "logit" for binomial).

- m.rate:

  A numeric vector of mismatch rates. If the length is 1, it is
  replicated for all blocks. If length \> 1, it must match the number of
  unique blocks.

- audit.size:

  a vector of block sizes in the audit sample (selected by simple random
  sampling) if used to estimate the m.rate (optional). If a single value
  is provided, assume the same value for all blocks and put out a
  warning.

- blocks:

  A vector indicating the block membership of each observation.

- weight.matrix:

  A character string specifying the weighting method ("ratio-type",
  "Lahiri-Larsen", "BLUE", or "all" (default))

- control:

  an optional list variable to of control arguments including
  "init.beta" for the initial outcome model coefficient estimates) - by
  default is the naive estimator when the weight matrix is ratio-type or
  Lahiri-Larsen and is the Lahiri-Larsen estimator for the BLUE weight
  matrix.

- ...:

  Pass control arguments directly.

## Value

A list of results:

- coefficients:

  A named vector (or matrix) of coefficients for the outcome model.

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

- rank:

  The numeric rank of the fitted linear model.

- df.residual:

  The residual degrees of freedom.

- df.null:

  The residual degrees of freedom for the null model.

- family:

  The `family` object used.

- call:

  The matched call.

## References

Chambers, R. (2009). Regression analysis of probability-linked data.
*Official Statistics Research Series*, 4, 1-15.

## Examples

``` r
data(brfss)
brfss <- na.omit(brfss)

x <- cbind(1, subset(brfss, select = c(Height,Physhlth,Menthlth,Exerany)))
y <- brfss$Weight

fit <- glmELE(x, y, family = "gaussian",
             m.rate = unique(brfss$m.rate), blocks = brfss$imonth,
             weight.matrix = "BLUE")
```
