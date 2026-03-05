# Fit a CoxPH Model Assuming Exchangeable Linkage Errors

Fit Cox proportional hazards regression adjusted for mismatched data
based on the approach developed in Vo et al., 2024 assuming exchangeable
linkage errors. Block-wise mismatch rates are assumed to be known.

## Usage

``` r
coxphELE(
  x,
  y,
  cens,
  m.rate,
  blocks,
  audit.size = NULL,
  control = list(init.beta = NULL),
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

- m.rate:

  block-wise mismatch rates (should be a vector with length equal to the
  number of blocks) - by default assume a single block.

- blocks:

  block indicators.

- audit.size:

  a vector of block sizes in the audit sample (selected by simple random
  sampling) if used to estimate the m.rate (optional). If a single value
  is provided, assume the same value for all blocks and put out a
  warning.

- control:

  an optional list variable to of control arguments including
  "init.beta" for the initial outcome model coefficient estimates) - by
  default is the naive estimator.

- ...:

  the option to directly pass "control" arguments

## Value

a list of results from the function called depending on the "family"
specified.

- coefficients:

  the outcome model coefficient estimates

- var:

  the variance-covariance matrix

- linear.predictors:

  the linear predictors

- means:

  Column means of the covariate matrix `x`.

- n:

  Number of observations.

- nevent:

  Number of events.

## References

Vo, T. H., Garès, V., Zhang, L. C., Happe, A., Oger, E., Paquelet, S., &
Chauvet, G. (2024). Cox regression with linked data. Statistics in
Medicine, 43(2), 296-314.  
