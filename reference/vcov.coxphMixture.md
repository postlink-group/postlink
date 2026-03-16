# Extract Variance-Covariance Matrix from a coxphMixture Object

Extracts the variance-covariance matrix of the main parameters from a
fitted `coxphMixture` object. The matrix is estimated using Louis'
method (1982) to account for the missing data structure (latent match
status) inherent in the mixture model.

## Usage

``` r
# S3 method for class 'coxphMixture'
vcov(object, ...)
```

## Arguments

- object:

  An object of class `coxphMixture`.

- ...:

  Additional arguments (currently ignored).

## Value

A matrix of the estimated covariances between the parameter estimates.
The rows and columns correspond to the outcome model coefficients
(`beta`) and the mismatch model coefficients (`gamma`).

## References

Louis, T. A. (1982). Finding the observed information matrix when using
the EM algorithm. *Journal of the Royal Statistical Society: Series B
(Methodological)*, 44(2), 226-233.

## Examples

``` r
library(survival)
set.seed(201)

# Simulate survival data (N = 200)
n <- 200
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.5)
true_time <- rexp(n, rate = exp(0.5 * x1 - 0.5 * x2))
cens_time <- rexp(n, rate = 0.5)

# Simulate auxiliary match scores and linkage errors
# Lower match scores correspond to a higher probability of mismatch
match_score <- runif(n, 0.5, 1.0)
is_mismatch <- rbinom(n, 1, prob = 1 - match_score)

# Induce linkage errors by shuffling covariates of mismatched records
linked_x1 <- x1
linked_x2 <- x2
mis_idx <- which(is_mismatch == 1)
shuffled_idx <- sample(mis_idx)
linked_x1[mis_idx] <- x1[shuffled_idx]
linked_x2[mis_idx] <- x2[shuffled_idx]

linked_data <- data.frame(
  time = pmin(true_time, cens_time),
  status = as.numeric(true_time <= cens_time),
  x1 = linked_x1, x2 = linked_x2,
  match_score = match_score
)

# Fit the Cox PH Mixture Model (Slawski et al., 2023)
adj <- adjMixture(linked.data = linked_data, m.formula = ~ match_score)
fit <- plcoxph(Surv(time, status) ~ x1 + x2, adjustment = adj,
               control = list(max.iter = 15))

# Extract the variance-covariance matrix
# Note: For outcome coefficients (beta) and mismatch coefficients (gamma)
vmat <- vcov(fit)
print(vmat)
#>                        x1          x2 m.(Intercept) m.match_score
#> x1             0.05642079 -0.04016553     0.1909437    -0.2567206
#> x2            -0.04016553  0.17945502    -0.2074337     0.2891386
#> m.(Intercept)  0.19094373 -0.20743369     5.4601921    -6.8925090
#> m.match_score -0.25672058  0.28913860    -6.8925090     8.8613058
```
