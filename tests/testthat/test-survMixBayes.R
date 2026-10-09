# Tests for the Bayesian survreg mixture engine (C++ Gibbs sampler)
#
# survregMixBayes() basic object validity on simple synthetic survival data for
# both supported distributions (gamma, weibull), with right censoring.

local_edition(3)

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

  # scale = exp(eta) (unit multiplicative scale parameter)
  scale_vec <- ifelse(z_true == 1, exp(eta1), exp(eta2))
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

expect_valid_surv_fit <- function(fit, X, dist) {
 expect_s3_class(fit, "survMixBayes")
 expect_equal(fit$dist, dist)

 expect_true(is.matrix(fit$m_samples))
 expect_equal(ncol(fit$m_samples), nrow(X))
 expect_true(all(fit$m_samples %in% c(1L, 2L)))

 expect_true(is.list(fit$estimates))
 expect_true(is.matrix(fit$estimates$coefficients))
 expect_true(is.matrix(fit$estimates$coefficients2))
 expect_equal(ncol(fit$estimates$coefficients), ncol(X))
 expect_equal(ncol(fit$estimates$coefficients2), ncol(X))
 expect_equal(colnames(fit$estimates$coefficients), colnames(X))
 expect_equal(colnames(fit$estimates$coefficients2), colnames(X))

 S <- nrow(fit$m_samples)
 expect_length(fit$estimates$theta, S)
 expect_true(all(fit$estimates$theta > 0 & fit$estimates$theta < 1))
 expect_length(fit$estimates$shape, S)
 expect_length(fit$estimates$shape2, S)
 expect_true(all(fit$estimates$shape > 0))
 expect_true(all(fit$estimates$shape2 > 0))
 if (dist == "weibull") {
  expect_length(fit$estimates$scale, S)
  expect_length(fit$estimates$scale2, S)
  expect_true(all(fit$estimates$scale > 0))
  expect_true(all(fit$estimates$scale2 > 0))
 } else {
  expect_null(fit$estimates$scale)
 }
}

# ------------------------------------------------------------------------------
# Test 1: gamma survival mixture
# ------------------------------------------------------------------------------
test_that("survregMixBayes runs MCMC and returns a valid object (gamma)", {
 skip_if_not_installed("survival")

 d <- generate_bayessurv_mixture_data(family = "gamma", seed = 321, N = 120)

 fit <- suppressMessages(survregMixBayes(
  X = d$X,
  y = d$y,
  dist = "gamma",
  control = list(
   iterations = 1500,
   burnin.iterations = 500,
   seed = 321
  )
 ))

 expect_valid_surv_fit(fit, d$X, "gamma")
 expect_equal(nrow(fit$m_samples), 1000L)

 # The slope of the correct-match component is well identified
 expect_equal(unname(colMeans(fit$estimates$coefficients))[2], d$true$beta1[2], tolerance = 0.3)
 expect_gt(mean(fit$estimates$theta), 0.4)
})

# ------------------------------------------------------------------------------
# Test 2: weibull survival mixture (matrix and list responses)
# ------------------------------------------------------------------------------
test_that("survregMixBayes runs MCMC and returns a valid object (weibull)", {
 d <- generate_bayessurv_mixture_data(family = "weibull", seed = 11, N = 120)
 y_mat <- cbind(time = d$dat$y, event = d$dat$status)

 fit <- suppressMessages(survregMixBayes(
  X = d$X,
  y = y_mat,
  dist = "Weibull",
  control = list(iterations = 1500, burnin.iterations = 500, seed = 5)
 ))
 expect_valid_surv_fit(fit, d$X, "weibull")
 expect_equal(unname(colMeans(fit$estimates$coefficients))[2], d$true$beta1[2], tolerance = 0.3)

 # list response gives the same result for the same seed
 fit2 <- suppressMessages(survregMixBayes(
  X = d$X,
  y = list(time = d$dat$y, event = d$dat$status),
  dist = "weibull",
  control = list(iterations = 1500, burnin.iterations = 500, seed = 5)
 ))
 expect_equal(fit2$estimates$coefficients, fit$estimates$coefficients)
})

# ------------------------------------------------------------------------------
# Test 3: safe matches and Path B
# ------------------------------------------------------------------------------
test_that("survregMixBayes supports safe matches and the logistic theta model", {
 d <- generate_bayessurv_mixture_data(family = "gamma", seed = 8, N = 120)
 safe <- as.integer(d$z_true == 1 & seq_len(120) %% 4 == 0)

 fit <- suppressMessages(survregMixBayes(d$X, d$y, dist = "gamma", safe.matches = safe,
                        control = list(iterations = 500, burnin.iterations = 200, seed = 2)))
 expect_valid_surv_fit(fit, d$X, "gamma")
 expect_true(all(fit$m_samples[, safe == 1L] == 1L))

 Z <- cbind(1, stats::rnorm(120)); colnames(Z) <- c("(Intercept)", "w")
 fitB <- suppressMessages(survregMixBayes(d$X, d$y, dist = "gamma", Z = Z,
                         control = list(iterations = 500, burnin.iterations = 200, seed = 2)))
 expect_s3_class(fitB, "survMixBayes")
 expect_true(fitB$use_logistic)
 expect_null(fitB$estimates$theta)
 expect_equal(dim(fitB$estimates$gamma), c(300L, 2L))
 expect_equal(colnames(fitB$estimates$gamma), colnames(Z))
})

# ------------------------------------------------------------------------------
# Test 4: input validation
# ------------------------------------------------------------------------------
test_that("survregMixBayes validates inputs", {
 d <- generate_bayessurv_mixture_data(family = "gamma", seed = 1, N = 40)
 ctrl <- list(iterations = 100, burnin.iterations = 20, seed = 1)
 expect_error(suppressMessages(survregMixBayes(d$X, d$y, dist = "lognormal", control = ctrl)), "gamma' or 'weibull")
 expect_error(suppressMessages(survregMixBayes(d$X, d$y[-1, ], dist = "gamma", control = ctrl)), "length nrow")
 expect_error(suppressMessages(survregMixBayes(d$X, d$y, dist = "gamma", safe.matches = 1L, control = ctrl)), "must have length")
})
