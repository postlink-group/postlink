# Tests for the Bayesian GLM mixture engine (C++ Gibbs sampler)
#
# glmMixBayes() basic object validity on simple synthetic mixture data for all
# supported families (gaussian, poisson, binomial, gamma), plus safe matches
# and the logistic match-probability model (Path B).

local_edition(3)

# ------------------------------------------------------------------------------
# Helper: generate simple synthetic 2-component mixture GLM data
# ------------------------------------------------------------------------------
generate_bayesglm_mixture_data <- function(
  family = "gaussian",
  seed = 123,
  n = 120,
  theta = 0.6
) {
 family <- tolower(trimws(as.character(family)[1]))
 if (!family %in% c("gaussian", "poisson", "binomial", "gamma")) {
  stop("family must be one of gaussian/poisson/binomial/gamma", call. = FALSE)
 }

 set.seed(seed)

 # Design: intercept + one continuous + one binary predictor
 x1 <- stats::runif(n, -2, 2)
 x2 <- stats::rbinom(n, 1, 0.5)
 X  <- cbind(1, x1, x2)
 colnames(X) <- c("(Intercept)", "x1", "x2")

 # Mixture membership: 1 (with probability theta) or 2
 z <- 2L - stats::rbinom(n, 1, theta)

 # Choose component-specific parameters with decent separation
 if (family == "gaussian") {
  beta1 <- c(0.2,  1.0, -0.6)
  beta2 <- c(1.5, -0.8,  0.9)
  sigma1 <- 0.7
  sigma2 <- 1.1

  eta1 <- drop(X %*% beta1)
  eta2 <- drop(X %*% beta2)
  mu   <- ifelse(z == 1L, eta1, eta2)
  sd   <- ifelse(z == 1L, sigma1, sigma2)
  y    <- stats::rnorm(n, mean = mu, sd = sd)

  truth <- list(beta1 = beta1, beta2 = beta2,
                dispersion1 = sigma1^2, dispersion2 = sigma2^2, theta = theta)

 } else if (family == "poisson") {
  beta1 <- c(-0.2, 0.4, -0.3)
  beta2 <- c( 0.8, 0.2,  0.5)

  eta1 <- drop(X %*% beta1)
  eta2 <- drop(X %*% beta2)
  lam  <- ifelse(z == 1L, exp(eta1), exp(eta2))
  y    <- stats::rpois(n, lambda = lam)

  truth <- list(beta1 = beta1, beta2 = beta2, theta = theta)

 } else if (family == "binomial") {
  beta1 <- c(-0.5, 0.8, -0.6)
  beta2 <- c( 0.7,-0.4,  0.9)

  eta1 <- drop(X %*% beta1)
  eta2 <- drop(X %*% beta2)
  p    <- ifelse(z == 1L, stats::plogis(eta1), stats::plogis(eta2))
  y    <- stats::rbinom(n, size = 1, prob = p)

  truth <- list(beta1 = beta1, beta2 = beta2, theta = theta)

 } else if (family == "gamma") {
  # Use mean = exp(eta), with component-specific shape (phi)
  beta1 <- c(0.1, 0.3, -0.2)
  beta2 <- c(0.6,-0.2,  0.4)
  phi1  <- 3.0
  phi2  <- 2.0

  eta1  <- drop(X %*% beta1)
  eta2  <- drop(X %*% beta2)
  mu    <- ifelse(z == 1L, exp(eta1), exp(eta2))

  shape <- ifelse(z == 1L, phi1, phi2)
  rate  <- shape / mu
  y     <- stats::rgamma(n, shape = shape, rate = rate)

  truth <- list(beta1 = beta1, beta2 = beta2,
                dispersion1 = 1/phi1, dispersion2 = 1/phi2, theta = theta)
 }

 list(X = X, y = y, z = z, truth = truth)
}

expect_valid_glm_fit <- function(fit, X, has_dispersion) {
 expect_s3_class(fit, "glmMixBayes")

 expect_true(is.matrix(fit$m_samples))
 expect_equal(ncol(fit$m_samples), nrow(X))
 expect_true(all(fit$m_samples %in% c(1L, 2L)))

 expect_true(is.list(fit$estimates))
 expect_true(is.matrix(fit$estimates$coefficients))
 expect_true(is.matrix(fit$estimates$coefficients2))
 expect_equal(ncol(fit$estimates$coefficients), ncol(X))
 expect_equal(ncol(fit$estimates$coefficients2), ncol(X))
 expect_equal(nrow(fit$estimates$coefficients), nrow(fit$m_samples))
 expect_equal(colnames(fit$estimates$coefficients), colnames(X))
 expect_equal(colnames(fit$estimates$coefficients2), colnames(X))
 expect_true(all(is.finite(fit$estimates$coefficients)))
 expect_true(all(is.finite(fit$estimates$coefficients2)))

 # mismatch-indicator model on the scale of glmMixture()'s m.coefficients
 expect_true(is.matrix(fit$estimates$m.coefficients))
 expect_equal(nrow(fit$estimates$m.coefficients), nrow(fit$m_samples))
 if (isTRUE(fit$use_logistic)) {
  expect_equal(fit$estimates$m.coefficients, -fit$estimates$gamma)
 } else {
  expect_equal(colnames(fit$estimates$m.coefficients), "(Intercept)")
  expect_equal(fit$estimates$m.coefficients[, 1], stats::qlogis(1 - fit$estimates$theta))
 }
 # per-record posterior probability of a correct match
 expect_length(fit$match.prob, nrow(X))
 expect_equal(unname(fit$match.prob), unname(colMeans(fit$m_samples == 1L)))

 expect_equal("dispersion" %in% names(fit$estimates), has_dispersion)
 expect_equal("dispersion2" %in% names(fit$estimates), has_dispersion)
 expect_false(any(c("m.dispersion", "m.shape", "m.scale") %in% names(fit$estimates)))
 if (has_dispersion) {
  expect_true(all(fit$estimates$dispersion > 0))
  expect_true(all(fit$estimates$dispersion2 > 0))
 }
 expect_true(is.list(fit$diagnostics))
 acc <- unlist(fit$diagnostics$accept)
 expect_true(all(is.na(acc) | (acc >= 0 & acc <= 1)))
}

# ------------------------------------------------------------------------------
# Test 1: gaussian fit recovers a well-separated mixture
# ------------------------------------------------------------------------------
test_that("glmMixBayes runs MCMC and returns a valid object (gaussian)", {
 expect_error(suppressMessages(glmMixBayes(X = matrix(1, 3, 1), y = rnorm(3), family = "badfamily")),
              "gaussian, poisson, binomial, or gamma")

 dat <- generate_bayesglm_mixture_data(family = "gaussian", seed = 123, n = 150, theta = 0.6)
 X <- dat$X
 y <- dat$y

 custom_priors <- list(
  beta1 = "normal(0,3)",
  beta2 = "normal(0,4)",
  sigma1 = "cauchy(0,1.5)",
  sigma2 = "cauchy(0,2)",
  theta = "beta(2,2)"
 )

 fit <- suppressMessages(glmMixBayes(
  X = X,
  y = y,
  family = "gaussian",
  priors = custom_priors,
  control = list(
   iterations = 3000,
   burnin.iterations = 1000,
   seed = 123
  )
 ))

 expect_valid_glm_fit(fit, X, has_dispersion = TRUE)
 expect_equal(nrow(fit$m_samples), 2000L)
 expect_length(fit$estimates$theta, 2000L)
 expect_true(all(fit$estimates$theta > 0 & fit$estimates$theta < 1))

 # Posterior means should be close to the oracle fits on the true components
 # (the components are well separated, so the mixture recovers the split)
 oracle1 <- unname(stats::coef(stats::lm(y[dat$z == 1L] ~ X[dat$z == 1L, -1])))
 oracle2 <- unname(stats::coef(stats::lm(y[dat$z == 2L] ~ X[dat$z == 2L, -1])))
 expect_equal(unname(colMeans(fit$estimates$coefficients)), oracle1, tolerance = 0.15)
 expect_equal(unname(colMeans(fit$estimates$coefficients2)), oracle2, tolerance = 0.15)
 expect_equal(mean(fit$estimates$theta), mean(dat$z == 1L), tolerance = 0.15)

 # Posterior allocations agree with the true membership for most records
 z_hat <- ifelse(colMeans(fit$m_samples == 1L) > 0.5, 1L, 2L)
 expect_gt(mean(z_hat == dat$z), 0.8)
})

# ------------------------------------------------------------------------------
# Test 2: all families produce valid objects
# ------------------------------------------------------------------------------
test_that("glmMixBayes supports poisson, binomial and gamma families", {
 for (fam in c("poisson", "binomial", "gamma")) {
  dat <- generate_bayesglm_mixture_data(family = fam, seed = 321, n = 120, theta = 0.6)
  # under the majority convention the draws are relabelled after sampling,
  # also when the component-specific priors differ (the binomial defaults):
  # on these weakly separated binomial data the chain visits both labellings,
  # which relabelling resolves without a warning
  fit <- suppressMessages(glmMixBayes(
   X = dat$X, y = dat$y, family = fam,
   control = list(iterations = 600, burnin.iterations = 200, seed = 1)
  ))
  expect_valid_glm_fit(fit, dat$X, has_dispersion = fam == "gamma")
  expect_equal(nrow(fit$m_samples), 400L)
  expect_equal(fit$diagnostics$orientation, "majority component")
  expect_true(fit$diagnostics$relabelled)
  expect_gte(mean(fit$m_samples == 1L), 0.5)
 }
})

# ------------------------------------------------------------------------------
# Test 3: reproducibility, thinning, family objects and control overrides
# ------------------------------------------------------------------------------
test_that("glmMixBayes is reproducible given a seed and honours thin / ...", {
 dat <- generate_bayesglm_mixture_data(family = "poisson", seed = 5, n = 80)

 f1 <- suppressMessages(glmMixBayes(dat$X, dat$y, family = stats::poisson(),
                   control = list(iterations = 400, burnin.iterations = 100, seed = 99)))
 f2 <- suppressMessages(glmMixBayes(dat$X, dat$y, family = "poisson",
                   control = list(iterations = 400, burnin.iterations = 100, seed = 99)))
 expect_equal(f1$estimates$coefficients, f2$estimates$coefficients)
 expect_equal(f1$m_samples, f2$m_samples)

 # dots override control; thinning reduces the number of stored draws
 f3 <- suppressMessages(glmMixBayes(dat$X, dat$y, family = "poisson",
                   control = list(iterations = 400, burnin.iterations = 100, seed = 99),
                   iterations = 500, thin = 4))
 expect_equal(nrow(f3$m_samples), 100L)

 # the caller's RNG stream is left untouched
 set.seed(2024); before <- runif(1)
 set.seed(2024)
 invisible(suppressMessages(glmMixBayes(dat$X, dat$y, family = "poisson",
                       control = list(iterations = 200, burnin.iterations = 50, seed = 1))))
 expect_equal(runif(1), before)
})

# ------------------------------------------------------------------------------
# Test 4: safe matches are fixed in component 1
# ------------------------------------------------------------------------------
test_that("glmMixBayes keeps safe matches in the correct-match component", {
 dat <- generate_bayesglm_mixture_data(family = "binomial", seed = 11, n = 120)
 safe <- as.integer(dat$z == 1L & seq_len(120) %% 3 == 0)
 expect_gt(sum(safe), 0)

 fit <- suppressMessages(glmMixBayes(dat$X, dat$y, family = "binomial", safe.matches = safe,
                    m.rate = 0.4, m.rate.sd = 0.1,
                    control = list(iterations = 600, burnin.iterations = 200, seed = 3)))
 expect_valid_glm_fit(fit, dat$X, has_dispersion = FALSE)
 expect_true(all(fit$m_samples[, safe == 1L] == 1L))
 expect_true(all(fit$match.prob[safe == 1L] == 1))
 expect_false(fit$use_logistic)
 expect_error(suppressMessages(glmMixBayes(dat$X, dat$y, family = "binomial", safe.matches = safe[-1])),
              "must have length")
 # the abbreviation `safe` still reaches the argument (partial matching)
 fit_abbrev <- suppressMessages(glmMixBayes(dat$X, dat$y, family = "binomial", safe = safe,
                    m.rate = 0.4, m.rate.sd = 0.1,
                    control = list(iterations = 600, burnin.iterations = 200, seed = 3)))
 expect_identical(fit_abbrev$m_samples, fit$m_samples)
})

# ------------------------------------------------------------------------------
# Test 5: logistic match-probability model (Path B)
# ------------------------------------------------------------------------------
test_that("glmMixBayes fits the logistic theta model when Z is supplied", {
 set.seed(77)
 n <- 200
 w <- stats::rnorm(n)
 z <- stats::rbinom(n, 1, stats::plogis(1 + 2 * w)) + 1L
 z <- 3L - z                         # z == 1 with probability inv_logit(1 + 2w)
 x1 <- stats::rnorm(n)
 X <- cbind(1, x1); colnames(X) <- c("(Intercept)", "x1")
 mu <- ifelse(z == 1L, 1 + 2 * x1, -2 - 2 * x1)
 y <- stats::rnorm(n, mu, 0.5)
 Z <- cbind(1, w); colnames(Z) <- c("(Intercept)", "w")

 fit <- suppressMessages(glmMixBayes(X, y, family = "gaussian", Z = Z,
                    control = list(iterations = 1500, burnin.iterations = 500, seed = 4)))

 expect_valid_glm_fit(fit, X, has_dispersion = TRUE)
 expect_true(fit$use_logistic)
 expect_null(fit$estimates$theta)
 expect_true(is.matrix(fit$estimates$gamma))
 expect_equal(dim(fit$estimates$gamma), c(1000L, 2L))
 expect_equal(colnames(fit$estimates$gamma), colnames(Z))
 expect_true(is.finite(fit$diagnostics$accept$gamma))

 # the slope on w should be recovered with the right sign and magnitude
 gm <- colMeans(fit$estimates$gamma)
 expect_gt(gm[["w"]], 0.8)
 expect_equal(unname(colMeans(fit$estimates$coefficients)), c(1, 2), tolerance = 0.2)

 expect_error(suppressMessages(glmMixBayes(X, y, family = "gaussian", Z = Z[-1, , drop = FALSE])),
              "same number of rows")
})

# ------------------------------------------------------------------------------
# Test 6: input validation
# ------------------------------------------------------------------------------
test_that("glmMixBayes validates inputs", {
 dat <- generate_bayesglm_mixture_data(family = "gaussian", seed = 1, n = 40)
 ctrl <- list(iterations = 100, burnin.iterations = 20, seed = 1)

 expect_error(suppressMessages(glmMixBayes(dat$X, dat$y[-1], family = "gaussian", control = ctrl)), "length nrow")
 expect_error(suppressMessages(glmMixBayes(dat$X, c(dat$y[-1], NA), family = "gaussian", control = ctrl)), "missing or infinite values")
 expect_error(suppressMessages(glmMixBayes(dat$X, rep(0:2, length.out = 40), family = "binomial", control = ctrl)), "binary")
 expect_error(suppressMessages(glmMixBayes(dat$X, abs(dat$y) + 0.5, family = "poisson", control = ctrl)), "coerced to integer")
 expect_error(suppressMessages(glmMixBayes(dat$X, dat$y, family = "gamma", control = ctrl)), "strictly positive")
 expect_error(suppressMessages(glmMixBayes(dat$X, dat$y, family = "gaussian",
                          control = list(iterations = 100, burnin.iterations = 100))),
              "burnin.iterations")
 expect_error(suppressMessages(glmMixBayes(dat$X, dat$y, family = "gaussian",
                          priors = list(sigma1 = "uniform(0, 10)"), control = ctrl)),
              "not supported")
})

# ------------------------------------------------------------------------------
# The sampler honours every distribution of the prior menu. Every record is a
# known correct match, so component 2 has no data (its posterior is the prior)
# and no relabelling happens.
# ------------------------------------------------------------------------------
test_that("the sampler honours every distribution of the prior menu", {
 set.seed(21)
 n <- 150
 x <- rnorm(n)
 X <- cbind("(Intercept)" = 1, x = x)
 y <- 1 + 2 * x + rnorm(n, sd = 0.7)
 safe <- rep(1L, n)
 ctl <- list(iterations = 1000, burnin.iterations = 300, seed = 7)
 # every record is safe by design: silence that (expected) warning only
 quiet_all_safe <- function(expr) withCallingHandlers(suppressMessages(expr), warning = function(w) {
  if (grepl("All records are flagged as safe matches", conditionMessage(w))) invokeRestart("muffleWarning")
 })
 expect_warning(suppressMessages(glmMixBayes(X, 1 + 2 * x, "gaussian", safe.matches = safe, control = ctl)),
                "All records are flagged as safe matches")
 for (pr in c("normal(2, 0.001)", "cauchy(2, 0.001)", "student_t(3, 2, 0.001)",
              "gamma(4e6, 2e6)", "lognormal(0.693147, 0.001)", "inv_gamma(1e6, 2e6)")) {
  f <- quiet_all_safe(glmMixBayes(X, y, "gaussian", priors = list(sigma2 = pr), safe.matches = safe, control = ctl))
  expect_lt(abs(median(sqrt(f$estimates$dispersion2)) - 2), 0.01, label = pr)
 }
 f <- quiet_all_safe(glmMixBayes(X, y, "gaussian", priors = list(sigma2 = "exponential(1000)"),
                                   safe.matches = safe, control = ctl))
 expect_lt(median(sqrt(f$estimates$dispersion2)), 0.005)

 tt <- rweibull(n, 2, exp(0.5 + 0.8 * x))
 for (pr in c("gamma(16e6, 4e6)", "lognormal(1.386294, 0.001)")) {
  f <- quiet_all_safe(survregMixBayes(X, cbind(tt, 1L), dist = "weibull", priors = list(scale2 = pr),
                                        safe.matches = safe, control = ctl))
  expect_lt(abs(median(f$estimates$scale2) - 4), 0.01, label = pr)
 }
 f <- quiet_all_safe(survregMixBayes(X, cbind(tt, 1L), dist = "weibull",
                                       priors = list(shape2 = "student_t(3, 3, 0.001)"), safe.matches = safe, control = ctl))
 expect_lt(abs(median(f$estimates$shape2) - 3), 0.01)

 yb <- rbinom(n, 1, plogis(0.5 + x))
 # the tight prior contradicts the data (slope about 1) and dominates the estimate
 expect_warning(f <- quiet_all_safe(glmMixBayes(X, yb, "binomial", priors = list(beta1 = "normal(3, 0.001)"),
                                                  safe.matches = safe, control = ctl)),
                "The supplied prior of 'x' in the correct-match component, normal(3, 0.001), conflicts with the data",
                fixed = TRUE)
 expect_lt(abs(mean(f$estimates$coefficients[, "x"]) - 3), 0.01)

 Zg <- cbind("(Intercept)" = 1, w = rnorm(n))
 f <- quiet_all_safe(glmMixBayes(X, y, "gaussian", Z = Zg, priors = list(gamma_slope = "normal(-2, 0.001)"),
                                   safe.matches = safe, control = ctl))
 expect_lt(abs(median(f$estimates$gamma[, "w"]) + 2), 0.01)

 f <- quiet_all_safe(glmMixBayes(X, y, "gaussian", m.rate = 0.3, m.rate.sd = 0.005, safe.matches = safe, control = ctl))
 expect_lt(abs(mean(f$estimates$theta) - 0.7), 0.005)
})
