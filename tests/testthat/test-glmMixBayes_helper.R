# tests/testthat/test-glmMixBayes_helpers.R
# Unit tests for internal Bayesian mixture helper utilities
# These tests are fast and do not run MCMC.

local_edition(3)

test_that("fill_defaults fills glm defaults and preserves user overrides", {
 res <- postlink:::fill_defaults(
  priors = list(beta1 = "normal(0,3)"),
  p_family = "gaussian",
  model_type = "glm"
 )

 expect_type(res, "list")
 expect_equal(res$beta1, "normal(0,3)")
 expect_equal(res$beta2, "normal(0,5)")
 expect_equal(res$sigma1, "cauchy(0,2.5)")
 expect_equal(res$sigma2, "cauchy(0,2.5)")
 expect_equal(res$theta, "beta(1,1)")
})

test_that("fill_defaults returns survival defaults and errors on invalid inputs", {
 res <- postlink:::fill_defaults(
  priors = NULL,
  p_family = "weibull",
  model_type = "survival"
 )

 expect_equal(res$beta1, "normal(0,2)")
 expect_equal(res$beta2, "normal(0,2)")
 expect_equal(res$shape1, "gamma(2,1)")
 expect_equal(res$shape2, "gamma(2,1)")
 expect_equal(res$scale1, "gamma(2,1)")
 expect_equal(res$scale2, "gamma(2,1)")
 expect_equal(res$theta, "beta(1,1)")

 expect_error(
  postlink:::fill_defaults(priors = 1, p_family = "gaussian", model_type = "glm"),
  "must be a named list"
 )

 expect_error(
  postlink:::fill_defaults(priors = list(), p_family = "badfamily", model_type = "glm"),
  "must be one of"
 )

 expect_error(
  postlink:::fill_defaults(priors = list(), p_family = "gaussian", model_type = "badtype"),
  "Unknown model_type"
 )
})

test_that("parse_prior_string parses valid priors and rejects malformed strings", {
 p1 <- postlink:::parse_prior_string("normal(0, 5)")
 expect_equal(p1$dist, "normal")
 expect_equal(p1$args, c(0, 5))

 p2 <- postlink:::parse_prior_string("beta(2,2)")
 expect_equal(p2$dist, "beta")
 expect_equal(p2$args, c(2, 2))

 expect_error(
  postlink:::parse_prior_string("normal 0,5"),
  "must be in format"
 )

 expect_error(
  postlink:::parse_prior_string("normal(a,5)"),
  "Could not parse numeric arguments"
 )
})

test_that("prepare_mixbayes_priors returns expected flat prior data for glm gaussian", {
 pri <- list(
  beta1 = "normal(1,2)",
  beta2 = "normal(3,4)",
  sigma1 = "cauchy(0,1.5)",
  sigma2 = "cauchy(0,2.5)",
  theta = "beta(2,3)"
 )

 out <- postlink:::prepare_mixbayes_priors(pri, family = "gaussian", model_type = "glm")

 expect_type(out, "list")
 expect_equal(out$prior_beta1_mu, 1)
 expect_equal(out$prior_beta1_sd, 2)
 expect_equal(out$prior_beta2_mu, 3)
 expect_equal(out$prior_beta2_sd, 4)
 expect_equal(out$prior_theta_alpha, 2)
 expect_equal(out$prior_theta_beta, 3)
 expect_equal(out$prior_sigma1_args, c(0, 1.5))
 expect_equal(out$prior_sigma2_args, c(0, 2.5))
 expect_equal(out$prior_sigma1_dist, "cauchy")
 expect_equal(out$prior_sigma2_dist, "cauchy")
})

test_that("prepare_mixbayes_priors handles survival gamma exponential priors correctly", {
 pri <- list(
  beta1 = "normal(0,5)",
  beta2 = "normal(0,5)",
  theta = "beta(1,1)",
  phi1 = "exponential(2)",
  phi2 = "exponential(3)"
 )

 out <- postlink:::prepare_mixbayes_priors(pri, family = "gamma", model_type = "survival")

 # exponential(rate) is converted to gamma(shape = 1, rate = rate)
 expect_equal(out$prior_phi1_args, c(1, 2))
 expect_equal(out$prior_phi2_args, c(1, 3))
 expect_equal(out$prior_phi1_dist, "gamma")
 expect_equal(out$prior_phi2_dist, "gamma")
})

test_that("prepare_mixbayes_priors errors when prior argument length is wrong", {
 pri <- list(
  beta1 = "normal(0,5)",
  beta2 = "normal(0,5)",
  theta = "beta(1,1)",
  sigma1 = "cauchy(0)",
  sigma2 = "cauchy(0,2.5)"
 )

 expect_error(
  postlink:::prepare_mixbayes_priors(pri, family = "gaussian", model_type = "glm"),
  "expected 2 arguments"
 )
})

test_that("build_engine_priors expands intercept/slope priors and encodes scalar priors", {
 flat <- postlink:::prepare_mixbayes_priors(
  list(intercept1 = "normal(1, 2)", beta1 = "normal(0, 3)",
       intercept2 = "normal(-1, 4)", beta2 = "normal(0, 0.5)",
       sigma1 = "cauchy(0, 2.5)", sigma2 = "lognormal(0, 1)", theta = "beta(2, 3)"),
  family = "gaussian", model_type = "glm"
 )
 eng <- postlink:::build_engine_priors(flat, "gaussian", "glm", K = 3L)

 expect_equal(eng$beta1_mean, c(1, 0, 0))
 expect_equal(eng$beta1_sd, c(2, 3, 3))
 expect_equal(eng$beta2_mean, c(-1, 0, 0))
 expect_equal(eng$beta2_sd, c(4, 0.5, 0.5))
 expect_equal(eng$theta, c(2, 3))
 expect_equal(eng$disp1, c(2, 0, 2.5, 0))   # cauchy -> code 2 (+ unused third argument)
 expect_equal(eng$disp2, c(5, 0, 1, 0))     # lognormal -> code 5
 expect_null(eng$gamma_mean)

 # Path B: gamma priors expanded over the M columns of Z
 flatB <- postlink:::prepare_mixbayes_priors(
  list(gamma_intercept = "normal(0.5, 1)", gamma_slope = "normal(0, 2)"),
  family = "poisson", model_type = "glm", use_logistic = TRUE
 )
 engB <- postlink:::build_engine_priors(flatB, "poisson", "glm", K = 2L, M = 3L)
 expect_equal(engB$gamma_mean, c(0.5, 0, 0))
 expect_equal(engB$gamma_sd, c(1, 2, 2))
 expect_null(engB$theta)
 expect_null(engB$disp1)

 # Weibull: shape and scale priors, exponential -> gamma(1, rate)
 flatW <- postlink:::prepare_mixbayes_priors(
  list(shape1 = "exponential(2)", shape2 = "student_t(3, 0, 2)", scale2 = "lognormal(0, 1)"),
  family = "weibull", model_type = "survival"
 )
 engW <- postlink:::build_engine_priors(flatW, "weibull", "survival", K = 2L)
 expect_equal(engW$disp1, c(3, 1, 2, 0))
 expect_equal(engW$disp2, c(6, 3, 0, 2))
 expect_equal(engW$scale2, c(5, 0, 1, 0))
 expect_equal(engW$scale1, c(3, 2, 1, 0))
 # the Weibull scale is sampled on the log scale jointly with beta and only
 # accepts log-concave priors
 expect_error(postlink:::prepare_mixbayes_priors(list(scale1 = "normal(1, 1)"),
                                                 family = "weibull", model_type = "survival"),
              "not supported")

 expect_error(postlink:::.scalar_prior_code("uniform"), "not supported")
 flat_bad <- flat; flat_bad$prior_beta1_sd <- -1
 expect_error(postlink:::build_engine_priors(flat_bad, "gaussian", "glm", K = 3L), "must be positive")
})

test_that("ecr_iterative_two aligns two-component allocations", {
 set.seed(1)
 S <- 60; N <- 40
 z_ref <- sample(1:2, N, replace = TRUE)
 z <- matrix(z_ref, S, N, byrow = TRUE)
 # sprinkle a few random disagreements, then swap the labels of some draws
 flip <- matrix(runif(S * N) < 0.05, S, N)
 z[flip] <- 3L - z[flip]
 swapped_draws <- seq(1, S, by = 3)
 z[swapped_draws, ] <- 3L - z[swapped_draws, ]

 swap <- postlink:::ecr_iterative_two(z)
 expect_type(swap, "logical")
 expect_length(swap, S)
 expect_true(all(swap[swapped_draws]))
 expect_false(any(swap[-swapped_draws]))

 # perfectly consistent allocations: nothing to swap
 expect_false(any(postlink:::ecr_iterative_two(matrix(z_ref, 10, N, byrow = TRUE))))
 expect_error(postlink:::ecr_iterative_two(1:5), "S x N matrix")
})

test_that("align_mixture_labels permutes every component-specific draw consistently", {
 set.seed(2)
 S <- 30; N <- 25
 z_ref <- rep(1:2, length.out = N)
 z <- matrix(z_ref, S, N, byrow = TRUE)
 swapped <- c(2, 5, 9)
 z[swapped, ] <- 3L - z[swapped, ]

 beta1 <- matrix(1, S, 2); beta2 <- matrix(2, S, 2)
 disp1 <- rep(10, S); disp2 <- rep(20, S)
 theta <- rep(0.7, S)

 out <- postlink:::align_mixture_labels(
  z, pairs = list(beta = list(beta1, beta2), disp = list(disp1, disp2)), theta = theta
 )
 expect_equal(which(out$swapped), swapped)
 expect_true(all(out$z == matrix(z_ref, S, N, byrow = TRUE)))
 expect_equal(out$pairs$beta[[1]][swapped, ], matrix(2, 3, 2))
 expect_equal(out$pairs$beta[[2]][swapped, ], matrix(1, 3, 2))
 expect_equal(out$pairs$beta[[1]][-swapped, ], matrix(1, S - 3, 2))
 expect_equal(out$pairs$disp[[1]][swapped], rep(20, 3))
 expect_equal(out$pairs$disp[[2]][swapped], rep(10, 3))
 expect_equal(out$theta[swapped], rep(0.3, 3))
 expect_equal(out$theta[-swapped], rep(0.7, S - 3))

 # gamma draws are negated on swapped iterations (Path B)
 gam <- matrix(c(1, -2), S, 2, byrow = TRUE)
 outB <- postlink:::align_mixture_labels(z, pairs = list(beta = list(beta1, beta2)), gamma = gam)
 expect_equal(outB$gamma[swapped, ], matrix(c(-1, 2), 3, 2, byrow = TRUE))
 expect_equal(outB$gamma[-swapped, ], matrix(c(1, -2), S - 3, 2, byrow = TRUE))

 # global swap when label 2 dominates
 z_dom <- 3L - z
 expect_message(outG <- postlink:::align_mixture_labels(z_dom, pairs = list(beta = list(beta1, beta2)), theta = theta,
                                                       verbose = TRUE),
                "Global label swap")
 expect_true(outG$flipped)
 expect_silent(postlink:::align_mixture_labels(z_dom, pairs = list(beta = list(beta1, beta2)), theta = theta))
 expect_true(all(outG$z == matrix(z_ref, S, N, byrow = TRUE)))
 expect_equal(outG$theta[-swapped], rep(0.3, S - 3))
})

test_that(".with_rng_seed seeds reproducibly and restores the caller's RNG state", {
 set.seed(10); a <- runif(1); rest_a <- runif(1)
 set.seed(10); a2 <- runif(1)
 x1 <- postlink:::.with_rng_seed(123, runif(1))
 expect_equal(runif(1), rest_a)          # stream continues where it left off
 x2 <- postlink:::.with_rng_seed(123, runif(1))
 expect_equal(x1, x2)
 expect_equal(a, a2)
 set.seed(123); expect_equal(x1, runif(1))
 expect_equal(postlink:::.with_rng_seed(NULL, 5), 5)
})

test_that(".mixbayes_control resolves dots > control > default", {
 expect_equal(postlink:::.mixbayes_control("thin", list(thin = 3), list(thin = 2), 1), 3)
 expect_equal(postlink:::.mixbayes_control("thin", list(), list(thin = 2), 1), 2)
 expect_equal(postlink:::.mixbayes_control("thin", list(), list(), 1), 1)
 expect_equal(postlink:::.mixbayes_control("thin", list(), NULL, 1), 1)
})

# --- Tests for new features: beta_to_logit_normal, use_logistic, m.rate ---

test_that("beta_to_logit_normal converts beta prior to logit-normal", {
 # beta(1,1) -> uniform on (0,1) -> mean = digamma(1) - digamma(1) = 0,
 #   sd = sqrt(2*trigamma(1))
 result <- postlink:::beta_to_logit_normal("beta(1, 1)")
 parsed <- postlink:::parse_prior_string(result)
 expect_equal(parsed$dist, "normal")
 expect_equal(parsed$args[1], 0, tolerance = 0.01)
 expect_true(parsed$args[2] > 1)

 # beta(70,30) -> mean approx digamma(70) - digamma(30)
 result2 <- postlink:::beta_to_logit_normal("beta(70, 30)")
 parsed2 <- postlink:::parse_prior_string(result2)
 expect_equal(parsed2$args[1], digamma(70) - digamma(30), tolerance = 0.01)
 expect_true(parsed2$args[2] > 0)
 expect_true(parsed2$args[2] < 1)

 # Error on non-beta distribution
 expect_error(postlink:::beta_to_logit_normal("normal(0,1)"),
              "expects a beta prior")
})

test_that("mrate_to_beta direct moment matching on probability scale", {
 # m.rate = 0.3, sigma_theta = 0.05
 out <- postlink:::mrate_to_beta(0.3, sigma_theta = 0.05)
 # E[theta] = a/(a+b) = 0.7
 expect_equal(out$alpha / (out$alpha + out$beta), 0.7, tolerance = 1e-8)
 # Var[theta] = ab / [(a+b)^2 (a+b+1)] = 0.05^2
 var_check <- out$alpha * out$beta /
              ((out$alpha + out$beta)^2 * (out$alpha + out$beta + 1))
 expect_equal(var_check, 0.05^2, tolerance = 1e-8)

 # Error if sigma_theta^2 >= m.rate * (1 - m.rate)
 expect_error(postlink:::mrate_to_beta(0.3, sigma_theta = 0.5),
              "must be strictly less than")
})

test_that("prepare_mixbayes_priors Path B derives the gamma intercept prior from m.rate", {
 # Path B with m.rate = 0.3, m.rate.sd = 0.05: the intercept prior is the exact
 # logit-scale moment match of the Beta prior that Path A would use
 out <- postlink:::prepare_mixbayes_priors(
  priors = NULL,
  family = "binomial",
  model_type = "glm",
  use_logistic = TRUE,
  m.rate = 0.3,
  m.rate.sd = 0.05
 )
 expect_true("prior_gamma_intercept_mu" %in% names(out))
 expect_true("prior_gamma_slope_mu" %in% names(out))
 expect_false("prior_theta_alpha" %in% names(out))

 bp <- postlink:::mrate_to_beta(0.3, 0.05)
 mom <- postlink:::logit_beta_moments(bp$alpha, bp$beta)
 expect_equal(out$prior_gamma_intercept_mu, unname(mom["mu"]), tolerance = 1e-8)
 expect_equal(out$prior_gamma_intercept_sd, unname(mom["sd"]), tolerance = 1e-8)
 expect_equal(out$prior_gamma_intercept_mu, qlogis(0.7), tolerance = 0.05)
 # the implied prior on theta has (approximately) the requested moments
 set.seed(1)
 th <- plogis(rnorm(2e5, out$prior_gamma_intercept_mu, out$prior_gamma_intercept_sd))
 expect_equal(mean(th), 0.7, tolerance = 0.01)
 expect_equal(sd(th), 0.05, tolerance = 0.05)
 # slopes default to normal(0, 2.5); gamma_slope / legacy gamma override it
 expect_equal(out$prior_gamma_slope_mu, 0)
 expect_equal(out$prior_gamma_slope_sd, 2.5)
 out2 <- postlink:::prepare_mixbayes_priors(list(gamma = "normal(0, 1)"), "binomial", "glm", use_logistic = TRUE)
 expect_equal(out2$prior_gamma_slope_sd, 1)
 expect_equal(out2$prior_gamma_intercept_sd, unname(postlink:::logit_beta_moments(1, 1)["sd"]))
})

test_that("prepare_mixbayes_priors Path A: m.rate -> Beta via direct moment matching", {
 out <- postlink:::prepare_mixbayes_priors(
  priors = NULL,
  family = "binomial",
  model_type = "glm",
  use_logistic = FALSE,
  m.rate = 0.3,
  m.rate.sd = 0.05
 )
 # E[theta] = alpha/(alpha+beta) = 0.7
 a <- out$prior_theta_alpha; b <- out$prior_theta_beta
 expect_equal(a / (a + b), 0.7, tolerance = 1e-8)
 expect_equal(a * b / ((a + b)^2 * (a + b + 1)), 0.05^2, tolerance = 1e-8)
})

test_that("prepare_mixbayes_priors: an explicit theta prior takes precedence over m.rate", {
 expect_message(out <- postlink:::prepare_mixbayes_priors(
  priors = list(theta = "beta(2, 3)"),
  family = "poisson",
  model_type = "glm",
  use_logistic = FALSE,
  m.rate = 0.3, m.rate.sd = 0.05
 ), "not used")
 expect_equal(out$prior_theta_alpha, 2)
 expect_equal(out$prior_theta_beta, 3)
 expect_false("prior_gamma_intercept_mu" %in% names(out))
 # even when the explicit prior equals the default string
 expect_message(out <- postlink:::prepare_mixbayes_priors(list(theta = "beta(1,1)"), "poisson", "glm",
                                                          m.rate = 0.3, m.rate.sd = 0.05), "not used")
 expect_equal(c(out$prior_theta_alpha, out$prior_theta_beta), c(1, 1))
})

test_that("prepare_mixbayes_priors validates prior strings and reports unknown entries", {
 expect_error(postlink:::prepare_mixbayes_priors(list(beta1 = "student_t(3,0,1)"), "gaussian", "glm"),
              "not supported")
 expect_error(postlink:::prepare_mixbayes_priors(list(beta1 = "normal(0,-1)"), "gaussian", "glm"),
              "must be positive")
 expect_error(postlink:::prepare_mixbayes_priors(list(theta = "normal(0,1)"), "gaussian", "glm"),
              "not supported")
 expect_error(postlink:::prepare_mixbayes_priors(NULL, "gaussian", "glm", m.rate = 1.2), "strictly between")
 expect_error(postlink:::prepare_mixbayes_priors(NULL, "gaussian", "glm", m.rate = 0.05, m.rate.sd = 0.3),
              "strictly less than")
 expect_warning(postlink:::prepare_mixbayes_priors(list(beta3 = "normal(0,1)"), "gaussian", "glm"),
                "unrecognised")
 expect_warning(postlink:::prepare_mixbayes_priors(list(gamma_slope = "normal(0,1)"), "gaussian", "glm"),
                "ignored")
 # exponential(rate) is the same prior as gamma(1, rate)
 out <- postlink:::prepare_mixbayes_priors(list(phi1 = "exponential(3)"), "gamma", "glm")
 expect_equal(out$prior_phi1_dist, "gamma")
 expect_equal(out$prior_phi1_args, c(1, 3))
})

test_that(".check_acceptance warns only for low acceptance rates", {
 expect_silent(postlink:::.check_acceptance(list(beta1 = 0.8, beta2 = 0.9, gamma = NA)))
 expect_warning(postlink:::.check_acceptance(list(beta1 = 0.05, beta2 = 0.9, gamma = NA)),
                "acceptance rate for the coefficients of component 1")
})

test_that("align_mixture_labels: known correct matches identify the labels (no relabelling)", {
 set.seed(3)
 S <- 30; N <- 12
 safe <- c(1L, 1L, rep(0L, N - 2L))
 # the sampler keeps safe records in component 1; 8 of 10 other records are mismatches
 z <- matrix(c(1L, 1L, rep(2L, 8), 1L, 1L), S, N, byrow = TRUE)
 beta1 <- matrix(1, S, 2); beta2 <- matrix(2, S, 2); theta <- rep(0.25, S)

 # label 2 dominates, but safe records identify the labels: the draws are
 # returned exactly as sampled
 expect_silent(out <- postlink:::align_mixture_labels(
  z, pairs = list(beta = list(beta1, beta2)), theta = theta, safe = safe))
 expect_identical(out$z, z)
 expect_identical(out$pairs$beta[[1L]], beta1)
 expect_equal(out$theta, theta)
 expect_false(any(out$swapped))
 expect_equal(out$orientation, "safe matches")
 expect_false(out$relabelled)
 expect_equal(out$ecr_share, 0)
 # relabelling is refused when the labels are identified
 expect_error(postlink:::align_mixture_labels(z, list(), safe = safe, relabel = TRUE), "majority convention")
 expect_error(postlink:::align_mixture_labels(z, list(), orientation = "match-rate prior", relabel = TRUE),
              "majority convention")

 # draws that look label-switched relative to the modal allocation are kept
 # (consistently with their parameter draws); the ECR share is computed on
 # the records not flagged as safe and warns only above 10% of the draws
 z2 <- z; z2[1:3, 3:N] <- 3L - z2[1:3, 3:N]          # 3 of 30 draws: 10%
 expect_silent(out2 <- postlink:::align_mixture_labels(
  z2, pairs = list(beta = list(beta1, beta2)), theta = theta, safe = safe))
 expect_equal(out2$ecr_share, 0.1)
 expect_identical(out2$z, z2)
 expect_false(any(out2$swapped))
 expect_true(all(out2$z[, 1:2] == 1L))
 z4 <- z; z4[1:4, 3:N] <- 3L - z4[1:4, 3:N]          # 4 of 30 draws: 13.3%
 expect_warning(out4 <- postlink:::align_mixture_labels(
  z4, pairs = list(beta = list(beta1, beta2)), theta = theta, safe = safe),
  "visited both labellings.*13.3%.*safe matches identify component 1")
 expect_identical(out4$z, z4)
 expect_equal(out4$ecr_share, 4 / 30)

 # every record a known correct match: no error, nothing to align
 z_all <- matrix(1L, S, N)
 out3 <- postlink:::align_mixture_labels(z_all, pairs = list(beta = list(beta1, beta2)),
                                         theta = theta, safe = rep(1L, N))
 expect_identical(out3$z, z_all)
 expect_false(any(out3$swapped))
 expect_equal(out3$ecr_share, 0)
 expect_true(is.na(out3$share1))
})

test_that("align_mixture_labels: orientation is decided after the ECR alignment", {
 # raw counts: label 2 in 40% of cells (no naive pre-swap), but after aligning
 # draw 3 to the modal allocation label 2 dominates -> everything is flipped
 z <- rbind(c(2L, 2L, 2L, 1L, 1L), c(2L, 2L, 2L, 1L, 1L), c(1L, 1L, 1L, 1L, 1L))
 beta1 <- matrix(seq_len(3), 3, 1); beta2 <- matrix(-seq_len(3), 3, 1); theta <- c(0.4, 0.4, 0.9)
 expect_message(out <- postlink:::align_mixture_labels(z, pairs = list(beta = list(beta1, beta2)), theta = theta,
                                                      verbose = TRUE),
                "Global label swap")
 expect_true(mean(out$z == 1L) > 0.5)
 expect_equal(out$theta, ifelse(out$swapped, 1 - theta, theta))
 expect_equal(out$pairs$beta[[1L]][out$swapped, ], beta2[out$swapped, ])
})
