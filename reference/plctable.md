# Analysis of Contingency Tables with Linkage Error Adjustment

`plctable` constructs a contingency table and adjusts the fitted model
for mismatch errors.

## Usage

``` r
plctable(
  formula,
  adjustment,
  subset,
  na.action,
  exclude = c(NA, NaN),
  control = list(),
  ...
)
```

## Arguments

- formula:

  a formula object with the left and right hand sides specifying the
  column and row variable of the flat table, respectively.

- adjustment:

  An object inheriting from class `"adjustment"`, or a `list` containing
  the necessary parameter specifications.

- subset:

  An optional vector specifying a subset of observations to be used.

- na.action:

  A function which indicates what should happen when the data contain
  NAs.

- exclude:

  Vector of values to be excluded when forming the table (passed to
  `xtabs`).

- control:

  A list of parameters for controlling the linkage error adjustment
  process.

- ...:

  Additional arguments passed to the internal function.

## Value

An object representing the fitted model. The specific class and
structure of the returned object depend directly on the `adjustment`
method provided:

- If `adjustment` is of class `adjMixture`, returns an object of class
  [`ctableMixture`](https://postlink-group.github.io/postlink/reference/ctableMixture.md).

## See also

[`adjMixture`](https://postlink-group.github.io/postlink/reference/adjMixture.md),
[`ctableMixture`](https://postlink-group.github.io/postlink/reference/ctableMixture.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(102)
n <- 400

# Simulate true categorical data
exposure <- sample(c("low", "high"), n, replace = TRUE)
# True relationship: high exposure -> higher chance of disease
prob_disease <- ifelse(exposure == "high", 0.7, 0.3)
true_outcome <- ifelse(runif(n) < prob_disease, "disease", "healthy")

# Induce linkage (mismatch) errors at a fixed overall rate
true_mismatch_rate <- 0.20
is_mismatch <- rbinom(n, 1, true_mismatch_rate)

obs_outcome <- true_outcome
mismatch_idx <- which(is_mismatch == 1)
# Shuffle outcomes for the mismatched records
obs_outcome[mismatch_idx] <- sample(obs_outcome[mismatch_idx])

linked_df <- data.frame(exposure, outcome = obs_outcome)

# Specify the Adjustment Method
adj <- adjMixture(
  linked.data = linked_df,
  m.rate = true_mismatch_rate
)

# Fit the adjusted contingency table model
fit <- plctable(
  ~ exposure + outcome,
  adjustment = adj
)
} # }
```
