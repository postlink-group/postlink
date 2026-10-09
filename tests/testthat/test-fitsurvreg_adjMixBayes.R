# tests/testthat/test-fitsurvreg_adjMixBayes.R
# Unit-style tests for internal dispatch + data alignment for Bayesian engine.
# This file DOES NOT run MCMC. It mocks survregMixBayes().

local_edition(3)

test_that("Basic Dispatch (Bayes): fitsurvreg dispatches to fitsurvreg.adjMixBayes and passes correct data to engine", {
 skip_if_not_installed("survival")

 # Capture arguments passed to the mocked engine
 mock_env <- new.env(parent = emptyenv())
 mock_env$args <- NULL

 # Mock survregMixBayes inside the package namespace to avoid real sampling
 testthat::local_mocked_bindings(
  survregMixBayes = function(X, y, dist, control, ...) {
   mock_env$args <- list(
    X = X, y = y, dist = dist, control = control, dots = list(...)
   )

   # Return a minimal valid object with class survMixBayes
   out <- list(
    m_samples = matrix(0, nrow = 2, ncol = nrow(X)),
    estimates = list(
     coefficients   = matrix(0, nrow = 2, ncol = ncol(X),
                             dimnames = list(NULL, colnames(X))),
     coefficients2  = matrix(0, nrow = 2, ncol = ncol(X),
                             dimnames = list(NULL, colnames(X)))
    )
   )
   class(out) <- "survMixBayes"
   out
  },
  .package = "postlink"
 )

 set.seed(11)
 n <- 10
 df <- data.frame(
  time  = rexp(n, rate = 0.2),
  event = rbinom(n, 1, 0.7),
  x1 = rnorm(n),
  x2 = rbinom(n, 1, 0.5)
 )
 rownames(df) <- paste0("id", seq_len(n))

 X <- stats::model.matrix(~ x1 + x2, df)
 rownames(X) <- rownames(df)

 y <- survival::Surv(df$time, df$event)

 # Minimal Bayes adjustment object
 adj <- postlink::adjMixBayes(linked.data = df)

 fit <- postlink:::fitsurvreg(
  x = X,
  y = y,
  dist = "gamma",
  adjustment = adj,
  control = list(iterations = 10, burnin.iterations = 5, seed = 11, cores = 1)
 )

 expect_s3_class(fit, "survMixBayes")
 expect_true(is.list(mock_env$args))

 # Engine should receive X and y
 expect_true(is.matrix(mock_env$args$X))
 expect_true(inherits(mock_env$args$y, "Surv"))
 expect_equal(nrow(mock_env$args$X), nrow(X))
 expect_equal(ncol(mock_env$args$X), ncol(X))

 # dist should pass through
 expect_true(is.character(mock_env$args$dist))
 expect_equal(mock_env$args$dist, "gamma")

 # control must be passed through
 expect_true(is.list(mock_env$args$control))
})

test_that("Row Alignment (Bayes): fitsurvreg.adjMixBayes subsets linked.data to match x rows", {
 skip_if_not_installed("survival")

 mock_env <- new.env(parent = emptyenv())
 mock_env$args <- NULL

 testthat::local_mocked_bindings(
  survregMixBayes = function(X, y, dist, control, ...) {
   mock_env$args <- list(X = X, y = y, dist = dist, control = control, dots = list(...))
   out <- list(
    m_samples = matrix(0, nrow = 2, ncol = nrow(X)),
    estimates = list(
     coefficients   = matrix(0, nrow = 2, ncol = ncol(X),
                             dimnames = list(NULL, colnames(X))),
     coefficients2  = matrix(0, nrow = 2, ncol = ncol(X),
                             dimnames = list(NULL, colnames(X)))
    )
   )
   class(out) <- "survMixBayes"
   out
  },
  .package = "postlink"
 )

 # Build linked.data with 6 rows, but X only uses a subset of 3 ids
 set.seed(2)
 df <- data.frame(
  time  = rexp(6, rate = 0.3),
  event = rbinom(6, 1, 0.6),
  x1 = rnorm(6),
  x2 = rbinom(6, 1, 0.5)
 )
 rownames(df) <- paste0("id", 1:6)

 keep_ids <- c("id1", "id3", "id6")
 df_sub <- df[keep_ids, , drop = FALSE]

 X <- stats::model.matrix(~ x1 + x2, df_sub)
 rownames(X) <- rownames(df_sub)

 y <- survival::Surv(df_sub$time, df_sub$event)

 adj <- postlink::adjMixBayes(linked.data = df)

 fit <- postlink:::fitsurvreg(
  x = X,
  y = y,
  dist = "gamma",
  adjustment = adj,
  control = list(iterations = 10, burnin.iterations = 5, seed = 2, cores = 1)
 )

 expect_s3_class(fit, "survMixBayes")
 expect_true(is.list(mock_env$args))

 # Most important: passed X corresponds to subset ordering in X
 expect_equal(rownames(mock_env$args$X), rownames(X))
})

test_that("Input Validation (Bayes): fitsurvreg errors if adjustment is not adjMixBayes", {
 skip_if_not_installed("survival")

 set.seed(3)
 n <- 5
 df <- data.frame(
  time  = rexp(n, rate = 0.4),
  event = rbinom(n, 1, 0.7),
  x1 = rnorm(n),
  x2 = rbinom(n, 1, 0.5)
 )
 X <- stats::model.matrix(~ x1 + x2, df)
 y <- survival::Surv(df$time, df$event)

 expect_error(
  postlink:::fitsurvreg(
   x = X,
   y = y,
   dist = "gamma",
   adjustment = list(not = "an adjustment"),
   control = list(iterations = 10, burnin.iterations = 5, seed = 3, cores = 1)
  ),
  "no applicable method"
 )
})

test_that("fitsurvreg.adjMixBayes passes Z, safe.matches, m.rate and m.rate.sd aligned with the modelled rows", {
 captured <- new.env(parent = emptyenv())
 testthat::local_mocked_bindings(
  survregMixBayes = function(X, y, dist, control, ...) {
   captured$args <- list(X = X, y = y, dots = list(...))
   structure(list(
    m_samples = matrix(1L, 2, nrow(X)),
    estimates = list(
     coefficients   = matrix(0, 2, ncol(X), dimnames = list(NULL, colnames(X))),
     coefficients2  = matrix(0, 2, ncol(X), dimnames = list(NULL, colnames(X)))
    ),
    dist = dist
   ), class = "survMixBayes")
  },
  .package = "postlink"
 )
 set.seed(2)
 n <- 30
 d <- data.frame(time = rexp(n) + 0.1, status = rbinom(n, 1, 0.7), x = rnorm(n), z = rnorm(n),
                 safe = rbinom(n, 1, 0.3) == 1)
 d$x[c(4, 11)] <- NA
 adj <- adjMixBayes(linked.data = d, m.formula = ~ z, safe.matches = safe, m.rate = 0.3, m.rate.sd = 0.05)
 plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull", adjustment = adj,
           control = list(iterations = 10, burnin.iterations = 5))
 keep <- !is.na(d$x)
 dots <- captured$args$dots
 expect_equal(as.vector(dots$Z), as.vector(cbind(1, d$z[keep])))
 expect_equal(dim(dots$Z), c(sum(keep), 2L))
 expect_identical(dots$safe.matches, as.integer(d$safe[keep]))
 expect_equal(dots$m.rate, 0.3)
 expect_equal(dots$m.rate.sd, 0.05)
 expect_equal(unname(captured$args$y[, "time"]), d$time[keep])
})
