# tests/testthat/test-survMixBayes_methods.R
# Methods tests for Bayesian survreg mixture (real Stan MCMC)
# Refactored into 4 test_that blocks + helper (glmMixture_methods style),
#
# NOTE:
# - Skipped on CRAN.
# - Opt-in by setting RUN_STAN_TESTS=true.

local_edition(3)

skip_if_no_stan <- function() {
 skip_on_cran()
 skip_if_not_installed("rstan")
 skip_if_not_installed("label.switching")
 skip_if_not_installed("survival")
 if (!identical(Sys.getenv("RUN_STAN_TESTS"), "true")) {
  skip("Set RUN_STAN_TESTS=true to run real Stan MCMC tests.")
 }
 rstan::rstan_options(auto_write = TRUE)
 options(mc.cores = 1L)
}

# -------------------------------------------------------------------------
# Helper: keep data generation consistent with your test-survMixBayes.R
# -------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Helper: simulate synthetic 2-component survival mixture with censoring
# ------------------------------------------------------------------------------
generate_bayessurv_mixture_data <- function(family = c("gamma", "weibull"),
                                            seed = 321,
                                            N = 120) {
 family <- match.arg(family)
 set.seed(seed)

 # Design: intercept + one covariate
 X <- cbind(1, stats::runif(N, -2, 2))
 colnames(X) <- c("X1", "X2")

 # True coefficients
 beta1 <- c(0.5, 1.2)
 beta2 <- c(1.5, -0.8)

 # True component labels: 60% in component 1
 z_true <- 1 + stats::rbinom(N, size = 1, prob = 0.4)  # P(z=1)=0.6, P(z=2)=0.4

 y <- numeric(N)

 if (family == "gamma") {
  phi1 <- 3.0
  phi2 <- 2.0

  eta1 <- as.vector(X %*% beta1)
  eta2 <- as.vector(X %*% beta2)
  mu1 <- exp(eta1)
  mu2 <- exp(eta2)

  shape_vec <- ifelse(z_true == 1, phi1, phi2)
  rate_vec  <- ifelse(z_true == 1, phi1 / mu1, phi2 / mu2)

  y <- stats::rgamma(N, shape = shape_vec, rate = rate_vec)

  true_par <- list(beta1 = beta1, beta2 = beta2, phi1 = phi1, phi2 = phi2)

 } else if (family == "weibull") {
  shape1 <- 2.0
  shape2 <- 1.5

  eta1 <- as.vector(X %*% beta1)
  eta2 <- as.vector(X %*% beta2)

  # Common choice: log(scale) depends on linear predictor
  scale_vec <- ifelse(z_true == 1, exp(-eta1), exp(-eta2))
  shape_vec <- ifelse(z_true == 1, shape1, shape2)

  y <- stats::rweibull(N, shape = shape_vec, scale = scale_vec)

  true_par <- list(beta1 = beta1, beta2 = beta2,
                   shape1 = shape1, shape2 = shape2,
                   scale1 = 1, scale2 = 1)
 }

 # Random right-censoring (keep censor rate moderate & stable)
 cmin <- as.numeric(stats::quantile(y, 0.4))
 cmax <- as.numeric(stats::quantile(y, 0.9))
 censoring_time <- stats::runif(N, min = cmin, max = cmax)

 status <- as.integer(y <= censoring_time)   # 1=event observed
 y_obs  <- pmin(y, censoring_time)

 dat <- data.frame(X1 = X[, 1], X2 = X[, 2], y = y_obs, status = status)

 list(
  dat = dat,
  X = X,
  y = survival::Surv(dat$y, dat$status),
  z_true = z_true,
  family = family,
  true = true_par
 )
}

# Fit ONCE and reuse across tests (avoid running Stan 4 times)
fit_once <- local({
 cache <- NULL
 function() {
  if (!is.null(cache)) return(cache)

  skip_if_no_stan()
  d <- generate_bayessurv_mixture_data(family = "gamma", seed = 321, N = 100)

  cache <<- survregMixBayes(
   X = d$X,
   y = d$y,
   dist = "gamma",
   control = list(
    iterations = 2000,
    burnin.iterations = 1000,
    seed = 101,
    cores = 1
   )
  )
  cache
 }
})

test_that("survMixBayes print() and summary() work on a real fitted object", {
 fit <- fit_once()

 out <- utils::capture.output(print(fit))
 expect_true(length(out) > 0)

 s <- summary(fit)
 expect_s3_class(s, "summary.survMixBayes")

 out2 <- utils::capture.output(print(s))
 expect_true(length(out2) > 0)
})

test_that("survMixBayes vcov() returns a p x p matrix", {
 fit <- fit_once()
 d <- generate_bayessurv_mixture_data(family = "gamma", seed = 321, N = 100)
 X <- d$X

 V <- stats::vcov(fit)
 expect_true(is.matrix(V))
 expect_equal(nrow(V), ncol(X))
 expect_equal(ncol(V), ncol(X))
})

test_that("survMixBayes confint() returns the expected list structure", {
 fit <- fit_once()
 d <- generate_bayessurv_mixture_data(family = "gamma", seed = 321, N = 100)
 X <- d$X

 CI <- stats::confint(fit)
 expect_true(is.list(CI))

 # coef blocks (matrices)
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

 # scalar block (numeric length-2)
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
})

test_that("survMixBayes predict() works and mi_with() methods exist + run", {
 fit <- fit_once()
 d <- generate_bayessurv_mixture_data(family = "gamma", seed = 321, N = 100)
 X <- d$X

 # predict(): accept newdata matrix and return list(component1, component2)
 newx <- X[1:5, , drop = FALSE]
 pr <- stats::predict(fit, newdata = newx)

 expect_true(is.list(pr))
 expect_true(all(c("component1", "component2") %in% names(pr)))
 expect_equal(length(pr$component1), nrow(newx))
 expect_equal(length(pr$component2), nrow(newx))

 # mi_with.survMixBayes exists
 wm <- getS3method("mi_with", "survMixBayes", optional = TRUE)
 expect_true(is.function(wm))

 pool <- mi_with(fit)

 expect_s3_class(pool, "mi_link_pool_survreg")
 expect_true(is.numeric(pool$p_component1))
 expect_equal(length(pool$p_component1), nrow(X))

 # posterior allocation probabilities should be in [0, 1]
 expect_true(all(pool$p_component1 >= 0 & pool$p_component1 <= 1, na.rm = TRUE))

 # pooled object should be printable
 out3 <- utils::capture.output(print(pool))
 expect_true(length(out3) > 0)
})
