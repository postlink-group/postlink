# Print Method for `adjELE` Objects

Provides a concise summary of the adjustment object created by
[`adjELE`](https://postlink-group.github.io/postlink/reference/adjELE.md),
including linkage error assumptions, blocking structure, and weight
estimation settings.

## Usage

``` r
# S3 method for class 'adjELE'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `adjELE`.

- digits:

  Integer; the number of significant digits to use when printing numeric
  values. Defaults to 3.

- ...:

  Additional arguments passed to methods.

## Value

Invisibly returns the input object `x`.

## Details

This method inspects the internal structure of the adjustment object. It
calculates summaries for mismatch rates and audit sizes (e.g.,
means/ranges) if they vary across blocks, providing a snapshot of the
error assumption complexity. It safely handles cases where the reference
data is missing or empty.

## Examples

``` r
data(brfss, package = "postlink")

adj_object <- adjELE(linked.data = brfss,
                    m.rate = unique(brfss$m.rate),
                    blocks = imonth,
                    weight.matrix = "BLUE")
print(adj_object)
#> 
#> --- Adjustment Object: Exchangeable Linkage Errors (Chambers, 2009) ---
#> 
#> * Linked Data:
#>     Observations:   2,000
#> * Specification:
#>     Weight Matrix:   BLUE
#>     Blocks:         12 distinct blocks
#>     Mismatch Rate:  Variable (Mean: 0.275, Range: 0.201 - 0.339)
#>     Audit Sample:   None (Using known rates)
#> 
```
