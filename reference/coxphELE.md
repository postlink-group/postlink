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

## Examples

``` r
set.seed(101)
n <- 250
# 1. Simulate true data
age <- rnorm(n, 60, 8)
treatment <- rbinom(n, 1, 0.5)
# Simulate true hazard
# Older age and lack of treatment increase hazard
true_hazard <- exp(0.05 * (age - 60) - 0.5 * treatment)
true_time <- rexp(n, rate = true_hazard)
cens_time <- rexp(n, rate = 0.1)
obs_time <- pmin(true_time, cens_time)
status <- as.numeric(true_time <= cens_time) # 1 = event, 0 = censored
# 2. Induce 15% exchangeable linkage errors
mis_idx <- sample(1:n, size = floor(0.15 * n))
linked_age <- age
linked_treatment <- treatment
if(length(mis_idx) > 0) {
  false_link_idx <- sample(1:n, size = length(mis_idx), replace = TRUE)
  linked_age[mis_idx] <- age[false_link_idx]
  linked_treatment[mis_idx] <- treatment[false_link_idx]
}

# 3. Prepare matrices for the internal routine
cens_ele <- 1 - status
X_mat <- cbind(age = linked_age, treatment = linked_treatment)

# 4. Fit the model directly
# cens = 1 for censored, 0 for event
fit_ele <- coxphELE(X = X_mat, y = obs_time, cens = cens_ele, m.rate = 0.15)
#> Error in coxphELE(X = X_mat, y = obs_time, cens = cens_ele, m.rate = 0.15): argument "x" is missing, with no default
print(fit_ele$coefficients)
#> Error: object 'fit_ele' not found
```
