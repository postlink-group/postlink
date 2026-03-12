# Secondary Analysis Constructor Based on Mixture Modeling

Specifies the linked data and information on the underlying record
linkage process for the "General Framework for Regression with
Mismatched Data" developed by Slawski et al. (2025). This framework uses
a mixture model for pairs of linked records whose two components reflect
distributions conditional on match status, i.e., correct match or
mismatch. Inference is based on composite likelihood and the EM
algorithm. Examples of information about the underlying record linkage
process that can be incorporated into the method if available are the
assumed overall mismatch rate, safe matches, predictors of match status,
or predicted probabilities of correct matches.

## Usage

``` r
adjMixture(
  linked.data = NULL,
  m.formula = ~1,
  m.rate = NULL,
  safe.matches = NULL
)
```

## Arguments

- linked.data:

  A data.frame containing the linked dataset. If `NULL`, the function
  attempts to resolve variables specified in `m.formula` from the
  environment.

- m.formula:

  A one-sided formula object for the mismatch indicator model, with the
  covariates on the right of "~". The default is an intercept-only model
  corresponding to a constant mismatch rate.

- m.rate:

  Numeric; an optional estimate (a proportion between 0 and 1) to
  constrain the global mismatch rate estimate. Defaults to `NULL`.

- safe.matches:

  A logical vector or an unquoted variable name found in `linked.data`;
  an indicator variable for safe matches (TRUE : record can be treated
  as a correct match and FALSE : record may be mismatched). The default
  is FALSE for all matches.

## Value

An object of class `c("adjMixture", "adjustment")`. To minimize memory
overhead, the underlying `linked.data` is stored by reference within an
environment inside this object.

## Details

The constructor assumes that any variables defined in `m.formula` and
`safe.matches` are in `linked.data` or in the same environment. Explicit
provision of `linked.data` is strongly recommended for reproducibility
and to ensure the adjustment object fully encapsulates the necessary
data for downstream model fitting.

## References

Slawski, M., West, B. T., Bukke, P., Wang, Z., Diao, G., & Ben-David, E.
(2025). A general framework for regression with mismatched data based on
mixture modelling. *Journal of the Royal Statistical Society Series A:
Statistics in Society*, 188(3), 896-919.
[doi:10.1093/jrsssa/qnae083](https://doi.org/10.1093/jrsssa/qnae083)

Bukke, P., Ben-David, E., Diao, G., Slawski, M., & West, B. T. (2025).
Cox Proportional Hazards Regression Using Linked Data: An Approach Based
on Mixture Modeling. In *Frontiers of Statistics and Data Science* (pp.
181-200). Singapore: Springer Nature Singapore.
[doi:10.1007/978-981-96-0742-6_8](https://doi.org/10.1007/978-981-96-0742-6_8)

Slawski, M., Diao, G., Ben-David, E. (2021). A pseudo-likelihood
approach to linear regression with partially shuffled data. *Journal of
Computational and Graphical Statistics*. 30(4), 991-1003.
[doi:10.1080/10618600.2020.1870482](https://doi.org/10.1080/10618600.2020.1870482)

## See also

- [`plglm()`](https://postlink-group.github.io/postlink/reference/plglm.md)
  for generalized linear regression modeling

- [`plcoxph()`](https://postlink-group.github.io/postlink/reference/plcoxph.md)
  for Cox proportional hazards regression modeling

- [`plctable()`](https://postlink-group.github.io/postlink/reference/plctable.md)
  for contingency table analysis

## Examples

``` r
# Load the LIFE-M demo dataset
data(lifem)

# Phase 1: Adjustment Specification
# We model the correct match indicator via logistic regression using
# name commonness scores (commf, comml) and a 5% expected mismatch rate.
adj_object <- adjMixture(
 linked.data = lifem,
 m.formula = ~ commf + comml,
 m.rate = 0.05,
 safe.matches = hndlnk
)
```
