# Methods tests for Bayesian survreg mixture (real Stan MCMC)

test_that("survMixBayes methods work on a real fitted object", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("label.switching")
  skip_if_not_installed("survival")

  if (!identical(Sys.getenv("RUN_STAN_TESTS"), "true")) {
    skip("Set RUN_STAN_TESTS=true to run real Stan MCMC tests.")
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
    dist = "gamma",
    control = list(iterations = 300, burnin.iterations = 150, seed = 101, cores = 1)
  )

  # print()
  out <- utils::capture.output(print(fit))
  expect_true(length(out) > 0)

  # summary()
  s <- summary(fit)
  expect_s3_class(s, "summary.survMixBayes")
  out2 <- utils::capture.output(print(s))
  expect_true(length(out2) > 0)

  # vcov()
  V <- stats::vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), ncol(X))
  expect_equal(ncol(V), ncol(X))

  # confint(): returns a list of credible intervals
  CI <- stats::confint(fit)
  expect_true(is.list(CI))

  # coef blocks (matrices)
  # (at least one of coef1/coef2 should exist in typical fits)
  expect_true(any(c("coef1", "coef2") %in% names(CI)))

  if ("coef1" %in% names(CI)) {
   expect_true(is.matrix(CI$coef1))
   expect_equal(nrow(CI$coef1), ncol(X))
   expect_equal(ncol(CI$coef1), 2L)
   expect_true(all(CI$coef1[, 1] <= CI$coef1[, 2], na.rm = TRUE))
  }

  if ("coef2" %in% names(CI)) {
   expect_true(is.matrix(CI$coef2))
   expect_equal(nrow(CI$coef2), ncol(X))
   expect_equal(ncol(CI$coef2), 2L)
   expect_true(all(CI$coef2[, 1] <= CI$coef2[, 2], na.rm = TRUE))
  }

  # scalar blocks (numeric length-2)
  expect_true("theta" %in% names(CI))
  expect_true(is.numeric(CI$theta))
  expect_equal(length(CI$theta), 2L)
  expect_true(CI$theta[1] <= CI$theta[2])

  for (nm in c("shape1", "shape2", "scale1", "scale2")) {
   if (nm %in% names(CI)) {
    expect_true(is.numeric(CI[[nm]]))
    expect_equal(length(CI[[nm]]), 2L)
    expect_true(CI[[nm]][1] <= CI[[nm]][2])
   }
  }

  # predict(): accept newx matrix (implementation-dependent output)
  newx <- X[1:5, , drop = FALSE]
  pr <- stats::predict(fit, newdata = newx)
  expect_true(is.list(pr))
  expect_true(all(c("component1", "component2") %in% names(pr)))
  expect_equal(length(pr$component1), nrow(newx))
  expect_equal(length(pr$component2), nrow(newx))

  # mi_with(): posterior allocation based pooling for survMixBayes
  mm <- getS3method("mi_with", "survMixBayes", optional = TRUE)
  expect_true(is.function(mm))

  pool <- mi_with(fit)

  expect_s3_class(pool, "mi_link_pool_survreg")
  expect_true(is.numeric(pool$p_component1))
  expect_equal(length(pool$p_component1), n)

  # posterior allocation probabilities should be in [0, 1]
  expect_true(all(pool$p_component1 >= 0 & pool$p_component1 <= 1, na.rm = TRUE))

  # pooled object should be printable
  out3 <- utils::capture.output(print(pool))
  expect_true(length(out3) > 0)
})
