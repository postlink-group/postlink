# Methods tests for Bayesian GLM mixture (real Stan MCMC)
# Mirrors the intent of test-glmMixture_methods.R but runs on glmMixBayes.

local_edition(3)

test_that("glmMixBayes methods work on a real fitted object", {
 skip_on_cran()
 skip_if_not_installed("rstan")
 skip_if_not_installed("label.switching")

 if (!identical(Sys.getenv("RUN_STAN_TESTS"), "true")) {
  skip("Set RUN_STAN_TESTS=true to run real Stan MCMC tests.")
 }

 rstan::rstan_options(auto_write = TRUE)
 options(mc.cores = 1L)

 set.seed(42)
 n <- 80
 dat <- data.frame(
  y  = rnorm(n),
  x1 = rnorm(n),
  x2 = rbinom(n, 1, 0.5)
 )

 X <- stats::model.matrix(y ~ x1 + x2, dat)
 yv <- dat$y

 fit <- glmMixBayes(
  X = X,
  y = yv,
  family = "gaussian",
  control = list(iterations = 300, burnin.iterations = 150, seed = 42, cores = 1)
 )

 # print()
 out <- utils::capture.output(print(fit))
 expect_true(length(out) > 0)

 # summary()
 s <- summary(fit)
 expect_s3_class(s, "summary.glmMixBayes")
 out2 <- utils::capture.output(print(s))
 expect_true(length(out2) > 0)

 # vcov()
 V <- stats::vcov(fit)
 expect_true(is.matrix(V))
 expect_equal(nrow(V), ncol(X))
 expect_equal(ncol(V), ncol(X))

 # confint()
 CI <- stats::confint(fit)
 expect_true(is.matrix(CI))
 expect_equal(nrow(CI), ncol(X))
 expect_equal(ncol(CI), 2L)
 expect_true(all(CI[, 1] <= CI[, 2]))

 # predict(): should accept newx as a matrix
 newx <- X[1:5, , drop = FALSE]
 pr <- stats::predict(fit, newx = newx, type = "link")

 # Allow either numeric vector OR list with $fit (depending on your predict impl)
 expect_true(is.numeric(pr) || is.list(pr))
 if (is.list(pr) && "fit" %in% names(pr)) {
  expect_equal(length(pr$fit), nrow(newx))
 }

 # ---------------------------
 # mi_with(): posterior allocation based pooling
 # ---------------------------

 # S3 method should exist
 mm <- getS3method("mi_with", "glmMixBayes", optional = TRUE)
 expect_true(is.function(mm))

 # run pooling (use the same dat used to build X/y)
 pool <- mi_with(
  fit,
  data = dat,
  formula = y ~ x1 + x2,
  quietly = TRUE
 )

 expect_s3_class(pool, "mi_link_pool")
 expect_s3_class(pool, "mi_link_pool_glm")
 expect_true(is.numeric(pool$coef))
 expect_true(is.matrix(pool$vcov))
 expect_equal(nrow(pool$vcov), length(pool$coef))
 expect_equal(ncol(pool$vcov), length(pool$coef))
})
