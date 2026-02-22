# Tests for Bayesian survreg mixture (real Stan MCMC)
# Mirrors the intent of test-survregMixture.R (if present) but runs survregMixBayes / fitsurvreg.adjMixBayes.

testthat::test_that("survregMixBayes runs real MCMC and returns a valid object", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")
  testthat::skip_if_not_installed("label.switching")
  testthat::skip_if_not_installed("survival")

  if (!identical(Sys.getenv("RUN_STAN_TESTS"), "true")) {
    testthat::skip("Set RUN_STAN_TESTS=true to run real Stan MCMC tests.")
  }

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 1L)

  set.seed(321)
  n <- 80
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)

  # Simulate log-normal survival with independent censoring
  linpred <- 0.2 + 0.6 * x1 - 0.4 * x2
  ttrue <- exp(linpred + rnorm(n, sd = 0.5))
  cens  <- rexp(n, rate = 0.25)

  time  <- pmin(ttrue, cens)
  event <- as.integer(ttrue <= cens)

  dat <- data.frame(time = time, event = event, x1 = x1, x2 = x2)

  X <- stats::model.matrix(~ x1 + x2, dat)
  y <- survival::Surv(dat$time, dat$event)

  fit <- survregMixBayes(
    X = X,
    y = y,
    dist = "lognormal",
    control = list(iterations = 400, burnin.iterations = 200, seed = 321, cores = 1)
  )

  testthat::expect_s3_class(fit, "survMixBayes")
  testthat::expect_true(is.matrix(fit$m_samples))
  testthat::expect_equal(ncol(fit$m_samples), nrow(X))

  testthat::expect_true(is.list(fit$estimates))
  testthat::expect_true(is.matrix(fit$estimates$coefficients))
  testthat::expect_true(is.matrix(fit$estimates$m.coefficients))
  testthat::expect_equal(ncol(fit$estimates$coefficients), ncol(X))
  testthat::expect_equal(ncol(fit$estimates$m.coefficients), ncol(X))
})

testthat::test_that("fitsurvreg dispatches to fitsurvreg.adjMixBayes and runs real MCMC", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")
  testthat::skip_if_not_installed("label.switching")
  testthat::skip_if_not_installed("survival")

  if (!identical(Sys.getenv("RUN_STAN_TESTS"), "true")) {
    testthat::skip("Set RUN_STAN_TESTS=true to run real Stan MCMC tests.")
  }

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 1L)

  set.seed(11)
  n <- 70
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)

  linpred <- -0.1 + 0.5 * x1 + 0.3 * x2
  ttrue <- exp(linpred + rnorm(n, sd = 0.6))
  cens  <- rexp(n, rate = 0.3)

  time  <- pmin(ttrue, cens)
  event <- as.integer(ttrue <= cens)

  dat <- data.frame(time = time, event = event, x1 = x1, x2 = x2)
  rownames(dat) <- paste0("id", seq_len(n))

  X <- stats::model.matrix(~ x1 + x2, dat)
  rownames(X) <- rownames(dat)

  y <- survival::Surv(dat$time, dat$event)

  adj <- adjMixBayes(linked.data = dat)

  fit <- fitsurvreg(
    x = X,
    y = y,
    dist = "lognormal",
    adjustment = adj,
    control = list(iterations = 300, burnin.iterations = 150, seed = 11, cores = 1)
  )

  testthat::expect_s3_class(fit, "survMixBayes")
  testthat::expect_true(is.matrix(fit$m_samples))
  testthat::expect_equal(ncol(fit$m_samples), nrow(X))
})
