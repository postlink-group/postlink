# tests/testthat/test-fitglm_adjMixBayes.R
# Unit-style tests for internal dispatch + data alignment for Bayesian engine.
# This file DOES NOT run MCMC. It mocks glmMixBayes().

local_edition(3)

test_that("Basic Dispatch (Bayes): fitglm dispatches to fitglm.adjMixBayes and passes correct data to engine", {

 # Capture arguments passed to the mocked engine
 mock_env <- new.env(parent = emptyenv())
 mock_env$args <- NULL

 # Mock glmMixBayes inside the package namespace to avoid real sampling
 testthat::local_mocked_bindings(
  glmMixBayes = function(X, y, family, control, ...) {
   mock_env$args <- list(
    X = X, y = y, family = family, control = control, dots = list(...)
   )

   # Return a minimal valid object with class glmMixBayes
   out <- list(
    m_samples = matrix(0, nrow = 2, ncol = nrow(X)),
    estimates = list(
     coefficients   = matrix(0, nrow = 2, ncol = ncol(X),
                             dimnames = list(NULL, colnames(X))),
     coefficients2  = matrix(0, nrow = 2, ncol = ncol(X),
                             dimnames = list(NULL, colnames(X))),
     dispersion     = rep(1, 2),
     dispersion2    = rep(1, 2)
    )
   )
   class(out) <- "glmMixBayes"
   out
  },
  .package = "postlink"
 )

 set.seed(1)
 n <- 10
 df <- data.frame(
  y  = rnorm(n),
  x1 = rnorm(n),
  x2 = rbinom(n, 1, 0.5)
 )
 rownames(df) <- paste0("id", seq_len(n))

 X <- stats::model.matrix(y ~ x1 + x2, df)
 yv <- df$y
 rownames(X) <- rownames(df)

 # Minimal Bayes adjustment object
 adj <- postlink::adjMixBayes(linked.data = df)

 fit <- postlink:::fitglm(
  x = X,
  y = yv,
  family = stats::gaussian(),
  adjustment = adj,
  control = list(iterations = 10, burnin.iterations = 5, seed = 1, cores = 1)
 )

 expect_s3_class(fit, "glmMixBayes")
 expect_true(is.list(mock_env$args))

 # Engine should receive X and y (possibly after internal checks)
 expect_true(is.matrix(mock_env$args$X))
 expect_true(is.numeric(mock_env$args$y))
 expect_equal(nrow(mock_env$args$X), nrow(X))
 expect_equal(ncol(mock_env$args$X), ncol(X))
 expect_equal(mock_env$args$y, yv)

 # family: depending on your implementation, it may pass "gaussian" or a family object
 expect_true(
  is.character(mock_env$args$family) ||
   inherits(mock_env$args$family, "family")
 )

 # control must be passed through
 expect_true(is.list(mock_env$args$control))
})

test_that("Row Alignment (Bayes): fitglm.adjMixBayes subsets linked.data to match x rows", {

 mock_env <- new.env(parent = emptyenv())
 mock_env$args <- NULL

 testthat::local_mocked_bindings(
  glmMixBayes = function(X, y, family, control, ...) {
   mock_env$args <- list(X = X, y = y, family = family, control = control, dots = list(...))
   out <- list(
    m_samples = matrix(0, nrow = 2, ncol = nrow(X)),
    estimates = list(
     coefficients   = matrix(0, nrow = 2, ncol = ncol(X),
                             dimnames = list(NULL, colnames(X))),
     coefficients2  = matrix(0, nrow = 2, ncol = ncol(X),
                             dimnames = list(NULL, colnames(X))),
     dispersion   = rep(1, 2),
     dispersion2  = rep(1, 2)
    )
   )
   class(out) <- "glmMixBayes"
   out
  },
  .package = "postlink"
 )

 # Build linked.data with 6 rows, but X only uses a subset of 3 ids
 set.seed(2)
 df <- data.frame(
  y  = rnorm(6),
  x1 = rnorm(6),
  x2 = rbinom(6, 1, 0.5)
 )
 rownames(df) <- paste0("id", 1:6)

 # Subset used for fitting
 keep_ids <- c("id1", "id3", "id6")
 df_sub <- df[keep_ids, , drop = FALSE]

 X <- stats::model.matrix(y ~ x1 + x2, df_sub)
 yv <- df_sub$y
 rownames(X) <- rownames(df_sub)

 adj <- postlink::adjMixBayes(linked.data = df)

 fit <- postlink:::fitglm(
  x = X,
  y = yv,
  family = stats::gaussian(),
  adjustment = adj,
  control = list(iterations = 10, burnin.iterations = 5, seed = 2, cores = 1)
 )

 expect_s3_class(fit, "glmMixBayes")
 expect_true(is.list(mock_env$args))

 # Most important: passed X/y correspond to the subset ordering in X
 expect_equal(rownames(mock_env$args$X), rownames(X))
 expect_equal(mock_env$args$y, yv)
})

test_that("Input Validation (Bayes): fitglm errors if adjustment is not adjMixBayes", {

 set.seed(3)
 df <- data.frame(
  y  = rnorm(5),
  x1 = rnorm(5),
  x2 = rbinom(5, 1, 0.5)
 )
 X <- stats::model.matrix(y ~ x1 + x2, df)
 yv <- df$y

 expect_error(
  postlink:::fitglm(
   x = X,
   y = yv,
   family = stats::gaussian(),
   adjustment = list(not = "an adjustment"),
   control = list(iterations = 10, burnin.iterations = 5, seed = 3, cores = 1)
  ),
  "no applicable method"
 )
})

test_that("fitglm.adjMixBayes passes Z, safe.matches, m.rate and m.rate.sd aligned with the modelled rows", {
 captured <- new.env(parent = emptyenv())
 testthat::local_mocked_bindings(
  glmMixBayes = function(X, y, family, control, ...) {
   captured$args <- list(X = X, y = y, dots = list(...))
   structure(list(
    m_samples = matrix(1L, 2, nrow(X)),
    estimates = list(
     coefficients   = matrix(0, 2, ncol(X), dimnames = list(NULL, colnames(X))),
     coefficients2  = matrix(0, 2, ncol(X), dimnames = list(NULL, colnames(X)))
    )
   ), class = "glmMixBayes")
  },
  .package = "postlink"
 )
 set.seed(1)
 n <- 30
 d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n), safe = rbinom(n, 1, 0.3) == 1)
 d$y[c(2, 9, 20)] <- NA
 adj <- adjMixBayes(linked.data = d, m.formula = ~ z, safe.matches = safe, m.rate = 0.2, m.rate.sd = 0.05)
 plglm(y ~ x, family = "gaussian", adjustment = adj, control = list(iterations = 10, burnin.iterations = 5))
 keep <- !is.na(d$y)
 dots <- captured$args$dots
 expect_equal(as.vector(dots$Z), as.vector(cbind(1, d$z[keep])))
 expect_equal(dim(dots$Z), c(sum(keep), 2L))
 expect_equal(colnames(dots$Z), c("(Intercept)", "z"))
 expect_identical(dots$safe.matches, as.integer(d$safe[keep]))
 expect_equal(dots$m.rate, 0.2)
 expect_equal(dots$m.rate.sd, 0.05)
 expect_equal(unname(captured$args$y), d$y[keep])
 expect_equal(rownames(captured$args$X), rownames(d)[keep])

 # without linkage information nothing is passed on
 plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(linked.data = d),
       control = list(iterations = 10, burnin.iterations = 5))
 dots <- captured$args$dots
 expect_null(dots$Z)
 expect_null(dots$safe.matches)
 expect_null(dots$m.rate)
 expect_null(dots$m.rate.sd)
})
