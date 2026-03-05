# Print Method for `adjMixture` Objects

Provides a concise summary of the adjustment object created by
[`adjMixture`](https://postlink-group.github.io/postlink/reference/adjMixture.md),
including dataset dimensions, model specifications, and constraints.

## Usage

``` r
# S3 method for class 'adjMixture'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `adjMixture`.

- digits:

  Integer; the number of significant digits to use when printing numeric
  values (e.g., mismatch rates). Defaults to 3.

- ...:

  Additional arguments passed to methods.

## Value

Invisibly returns the input object `x`.

## Details

This method inspects the reference-based data environment to report the
number of linked records without copying the full dataset. It considers
cases where components (like constraints or safe matches) are
unspecified.

## Examples

``` r
# Load the LIFE-M demo dataset
data(lifem)

# Phase 1: Adjustment Specification
adj_object <- adjMixture(
 linked.data = lifem,
 m.formula = ~ commf + comml,
 m.rate = 0.05,
 safe.matches = hndlnk
)

# Print specified adjustment
print(adj_object)
#> 
#> --- Adjustment Object: Mixture Model (Slawski et al., 2025) ---
#> 
#> * Linked Data:
#>     Observations: 3,238
#> * Specification:
#>     Mismatch Model:        ~commf + comml
#>     Global Mismatch Rate: 0.05 (User Constrained)
#>     Safe Matches:         2,159 (66.7%)
#> 
```
