# Secondary Analysis Constructor Assuming ELE

Specifies the linked data and information on the underlying record
linkage process for regression of linked data assuming exchangeable
linkage errors (ELE) as developed by Chambers (2009) and Vo et al.
(2024). These approaches correct for bias from mismatch error via
weighting matrices estimated using known mismatch rates or clerical
reviews (audit samples).

## Usage

``` r
adjELE(
  linked.data,
  m.rate,
  audit.size = NULL,
  blocks = NULL,
  weight.matrix = c("ratio", "LL", "BLUE", "all")
)
```

## Arguments

- linked.data:

  A data.frame containing the data after record linkage.

- m.rate:

  Numeric vector; known or estimated probability of mismatch for each
  record or block. Values must be between 0 and 1. Can be a single
  global rate, a vector of length equal to the number of unique blocks,
  or a vector of length equal to the number of rows in `linked.data`.

- audit.size:

  Numeric vector; If the m.rate is estimated, provide sample sizes for
  the clerical review audit. Used for variance estimation. If provided,
  must align with `blocks` similar to `m.rate`. Defaults to `NULL`
  (assume m.rate is known).

- blocks:

  A vector or an unquoted variable name found in `linked.data`
  identifying the blocking structure used during linkage. If `NULL`
  (default), all records are assumed to belong to a single block.

- weight.matrix:

  Character; the method for estimating the weight matrix. Must be one of
  "ratio" (default), "LL", "BLUE", or "all".

## Value

An object of class `c("adjELE", "adjustment")`. To minimize memory
overhead, the underlying `linked.data` is stored by reference within an
environment inside this object.

## Details

The constructor validates consistency between the mismatch rates, audit
sizes, and block identifiers. If `blocks` are provided, `m.rate` must be
specified either per-block (length equals number of unique blocks) or
per-record (length equals number of rows).

Explicit provision of `linked.data` is strongly recommended for
reproducibility and to ensure the adjustment object fully encapsulates
the necessary data for downstream model fitting.

## References

Chambers, R. (2009). Regression analysis of probability-linked data.
*Official Statistics Research Series*, 4, 1-15.

Vo, T. H., Garès, V., Zhang, L. C., Happe, A., Oger, E., Paquelet, S., &
Chauvet, G. (2024). Cox regression with linked data. *Statistics in
Medicine*, 43(2), 296-314.
[doi:10.1002/sim.9960](https://doi.org/10.1002/sim.9960)

## See also

- [`plglm()`](https://postlink-group.github.io/postlink/reference/plglm.md)
  for generalized linear regression modeling

- [`plcoxph()`](https://postlink-group.github.io/postlink/reference/plcoxph.md)
  for Cox proportional hazards regression modeling

## Examples

``` r
data(brfss, package = "postlink")

adj_object <- adjELE(linked.data = brfss,
                    m.rate = unique(brfss$m.rate),
                    blocks = imonth,
                    weight.matrix = "BLUE")
```
