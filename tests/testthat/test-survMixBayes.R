# Tests for Bayesian survreg mixture (real Stan MCMC)
# Mirrors the intent of test-survregMixture.R (if present) but runs survregMixBayes / fitsurvreg.adjMixBayes.

test_that("survregMixBayes runs real MCMC and returns a valid object", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("label.switching")
  skip_if_not_installed("survival")

  if (!identical(Sys.getenv("RUN_STAN_TESTS"), "true")) {
    skip("Set RUN_STAN_TESTS=true to run real Stan MCMC tests.")
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
    dist = "gamma",
    control = list(iterations = 400, burnin.iterations = 200, seed = 321, cores = 1)
  )

  expect_s3_class(fit, "survMixBayes")
  expect_true(is.matrix(fit$m_samples))
  expect_equal(ncol(fit$m_samples), nrow(X))

  expect_true(is.list(fit$estimates))
  expect_true(is.matrix(fit$estimates$coefficients))
  expect_true(is.matrix(fit$estimates$m.coefficients))
  expect_equal(ncol(fit$estimates$coefficients), ncol(X))
  expect_equal(ncol(fit$estimates$m.coefficients), ncol(X))
})

test_that("fitsurvreg dispatches to fitsurvreg.adjMixBayes and runs real MCMC", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("label.switching")
  skip_if_not_installed("survival")

  if (!identical(Sys.getenv("RUN_STAN_TESTS"), "true")) {
    skip("Set RUN_STAN_TESTS=true to run real Stan MCMC tests.")
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
    dist = "gamma",
    adjustment = adj,
    control = list(iterations = 300, burnin.iterations = 150, seed = 11, cores = 1)
  )

  expect_s3_class(fit, "survMixBayes")
  expect_true(is.matrix(fit$m_samples))
  expect_equal(ncol(fit$m_samples), nrow(X))
})
