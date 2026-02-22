# Methods tests for Bayesian survreg mixture (real Stan MCMC)

testthat::test_that("survMixBayes methods work on a real fitted object", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")
  testthat::skip_if_not_installed("label.switching")
  testthat::skip_if_not_installed("survival")

  if (!identical(Sys.getenv("RUN_STAN_TESTS"), "true")) {
    testthat::skip("Set RUN_STAN_TESTS=true to run real Stan MCMC tests.")
  }

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 1L)

  set.seed(101)
  n <- 80
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)

  linpred <- 0.1 + 0.4 * x1 - 0.2 * x2
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
    control = list(iterations = 300, burnin.iterations = 150, seed = 101, cores = 1)
  )

  # print()
  out <- utils::capture.output(print(fit))
  testthat::expect_true(length(out) > 0)

  # summary()
  s <- summary(fit)
  testthat::expect_s3_class(s, "summary.survMixBayes")
  out2 <- utils::capture.output(print(s))
  testthat::expect_true(length(out2) > 0)

  # vcov()
  V <- stats::vcov(fit)
  testthat::expect_true(is.matrix(V))
  testthat::expect_equal(nrow(V), ncol(X))
  testthat::expect_equal(ncol(V), ncol(X))

  # confint()
  CI <- stats::confint(fit)
  testthat::expect_true(is.matrix(CI))
  testthat::expect_equal(nrow(CI), ncol(X))
  testthat::expect_equal(ncol(CI), 2L)
  testthat::expect_true(all(CI[, 1] <= CI[, 2]))

  # predict(): accept newx matrix (implementation-dependent output)
  newx <- X[1:5, , drop = FALSE]
  pr <- stats::predict(fit, newx = newx, type = "lp")
  testthat::expect_true(is.numeric(pr) || is.list(pr))
  if (is.list(pr) && "fit" %in% names(pr)) {
    testthat::expect_equal(length(pr$fit), nrow(newx))
  }
})
