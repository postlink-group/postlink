# Tests for Bayesian GLM mixture (real Stan MCMC)
# Mirrors the intent of test-glmMixture.R but runs glmMixBayes / fitglm.adjMixBayes.
#
# NOTE:
# - These tests are skipped on CRAN.
# - Opt-in by setting environment variable RUN_STAN_TESTS=true.
#
# Example:
#   Sys.setenv(RUN_STAN_TESTS = "true")
#   devtools::test()

testthat::test_that("glmMixBayes runs real MCMC and returns a valid object", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")
  testthat::skip_if_not_installed("label.switching")

  if (!identical(Sys.getenv("RUN_STAN_TESTS"), "true")) {
    testthat::skip("Set RUN_STAN_TESTS=true to run real Stan MCMC tests.")
  }

  # Make Stan compilation more pleasant in local/CI runs
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 1L)

  set.seed(123)
  n <- 80
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  y  <- 0.5 + 1.2 * x1 - 0.8 * x2 + rnorm(n, sd = 0.7)

  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  X <- stats::model.matrix(y ~ x1 + x2, dat)
  yv <- dat$y

  fit <- glmMixBayes(
    X = X,
    y = yv,
    family = "gaussian",
    control = list(iterations = 400, burnin.iterations = 200, seed = 123, cores = 1)
  )

  testthat::expect_s3_class(fit, "glmMixBayes")
  testthat::expect_true(is.matrix(fit$m_samples))
  testthat::expect_equal(ncol(fit$m_samples), nrow(X))

  testthat::expect_true(is.list(fit$estimates))
  testthat::expect_true(is.matrix(fit$estimates$coefficients))
  testthat::expect_true(is.matrix(fit$estimates$m.coefficients))
  testthat::expect_equal(ncol(fit$estimates$coefficients), ncol(X))
  testthat::expect_equal(ncol(fit$estimates$m.coefficients), ncol(X))

  # Should carry coefficient names from X if available
  if (!is.null(colnames(X))) {
    testthat::expect_equal(colnames(fit$estimates$coefficients), colnames(X))
    testthat::expect_equal(colnames(fit$estimates$m.coefficients), colnames(X))
  }

  # Dispersion should be present for gaussian (and gamma)
  testthat::expect_true("dispersion" %in% names(fit$estimates))
  testthat::expect_true("m.dispersion" %in% names(fit$estimates))
})

testthat::test_that("fitglm dispatches to fitglm.adjMixBayes and runs real MCMC", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")
  testthat::skip_if_not_installed("label.switching")

  if (!identical(Sys.getenv("RUN_STAN_TESTS"), "true")) {
    testthat::skip("Set RUN_STAN_TESTS=true to run real Stan MCMC tests.")
  }

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 1L)

  set.seed(1)
  n <- 70
  dat <- data.frame(
    y  = rnorm(n),
    x1 = rnorm(n),
    x2 = rbinom(n, 1, 0.5)
  )

  # Stable rownames are important for the x <-> linked.data matching logic
  rownames(dat) <- paste0("id", seq_len(n))

  X <- stats::model.matrix(y ~ x1 + x2, dat)
  yv <- dat$y
  rownames(X) <- rownames(dat)

  adj <- adjMixBayes(linked.data = dat)

  fit <- fitglm(
    x = X,
    y = yv,
    family = stats::gaussian(),
    adjustment = adj,
    control = list(iterations = 300, burnin.iterations = 150, seed = 99, cores = 1)
  )

  testthat::expect_s3_class(fit, "glmMixBayes")
  testthat::expect_true(is.matrix(fit$m_samples))
  testthat::expect_equal(ncol(fit$m_samples), nrow(X))
})
