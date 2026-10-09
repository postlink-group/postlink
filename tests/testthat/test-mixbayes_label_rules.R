# Label-switching rules of the Bayesian mixture models (section 'Label
# switching' of ?glmMixBayes): what identifies component 1 is decided before
# sampling (safe matches > an asymmetric prior on the match probability > the
# majority convention). Under the match-rate prior the sampler orients its
# starting state towards the more probable labelling and a warning reports a
# weak identification; under the majority convention it orients its starting
# state so that component 1 is the majority (unless the component-specific
# priors make that labelling much less probable, which check chains confirm),
# and the draws are relabelled after sampling unless that exchange was
# refused, also when the component-specific priors differ.
local_edition(3)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Collect the warnings of an expression (muffled) together with its value.
with_warnings <- function(expr) {
 ws <- character()
 v <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
  ws <<- c(ws, conditionMessage(w))
  invokeRestart("muffleWarning")
 })
 list(value = v, warnings = ws)
}

# Log posterior of the reported draws of a gaussian or binomial glmMixBayes()
# fit without linkage covariates or safe matches, computed from the priors
# stored in the fit as the sampler computes diagnostics$lp (likelihood with the
# indicators integrated out, unnormalised priors on the natural scale). It
# equals diagnostics$lp only if the reported draws were sampled under these
# priors with these labels (or the labels were exchanged under exchangeable
# priors, which leaves the log posterior unchanged; under priors that differ,
# diagnostics$lp is that of the relabelled draw).
lp_reported <- function(f, X, y) {
 pr <- f$priors
 args <- function(key) postlink:::parse_prior_string(pr[[key]])$args
 K <- ncol(X)
 mean_sd <- function(i, b) list(m = c(i[1], rep(b[1], K - 1)), s = c(i[2], rep(b[2], K - 1)))
 p1 <- mean_sd(args("intercept1"), args("beta1"))
 p2 <- mean_sd(args("intercept2"), args("beta2"))
 ab <- args("theta")
 B1 <- f$estimates$coefficients; B2 <- f$estimates$coefficients2; th <- f$estimates$theta
 vapply(seq_len(nrow(B1)), function(s) {
  e1 <- drop(X %*% B1[s, ]); e2 <- drop(X %*% B2[s, ])
  if (f$family == "gaussian") {
   sd1 <- sqrt(f$estimates$dispersion[s]); sd2 <- sqrt(f$estimates$dispersion2[s])
   l1 <- stats::dnorm(y, e1, sd1, log = TRUE); l2 <- stats::dnorm(y, e2, sd2, log = TRUE)
   c1 <- args("sigma1"); c2 <- args("sigma2")          # half-cauchy(loc, scale)
   lpd <- -log1p(((sd1 - c1[1]) / c1[2])^2) - log1p(((sd2 - c2[1]) / c2[2])^2)
  } else {
   l1 <- stats::dbinom(y, 1, stats::plogis(e1), log = TRUE)
   l2 <- stats::dbinom(y, 1, stats::plogis(e2), log = TRUE)
   lpd <- 0
  }
  a <- log(th[s]) + l1; b <- log1p(-th[s]) + l2
  sum(pmax(a, b) + log1p(exp(-abs(a - b)))) -
   0.5 * sum(((B1[s, ] - p1$m) / p1$s)^2) - 0.5 * sum(((B2[s, ] - p2$m) / p2$s)^2) +
   (ab[1] - 1) * log(th[s]) + (ab[2] - 1) * log1p(-th[s]) + lpd
 }, numeric(1))
}

# Engine priors of a model, from user priors.
engine_priors <- function(priors, family, model_type = "glm", use_logistic = FALSE, K = 2L, M = 0L,
                          m.rate = NULL) {
 pf <- postlink:::prepare_mixbayes_priors(priors, family, model_type, use_logistic = use_logistic,
                                          m.rate = m.rate)
 postlink:::build_engine_priors(pf, family, model_type, K = K, M = M)
}

# ------------------------------------------------------------------------------
# The rule is decided from the engine priors and the safe matches
# ------------------------------------------------------------------------------

test_that(".label_rule() decides the orientation from the safe matches and the engine priors", {
 rule <- function(...) postlink:::.label_rule(engine_priors(...), integer(10))
 # exchangeable defaults: majority convention, oriented before and relabelled after sampling
 r <- rule(NULL, "gaussian")
 expect_equal(r$orientation, "majority component")
 expect_equal(r$pre_orient, "majority")
 expect_true(r$exchangeable)
 for (fam in c("poisson", "gamma")) expect_true(rule(NULL, fam)$exchangeable)
 for (dist in c("gamma", "weibull")) expect_true(rule(NULL, dist, "survival")$exchangeable)
 # binomial defaults differ between the components (beta1 normal(0, 2.5), beta2 normal(0, 5))
 r <- rule(NULL, "binomial")
 expect_equal(r$orientation, "majority component")
 expect_equal(r$pre_orient, "majority")
 expect_false(r$component_priors_identical)
 expect_false(r$exchangeable)
 # any component-specific prior that differs makes the components non-exchangeable
 expect_false(rule(list(intercept2 = "normal(0, 9)"), "gaussian")$exchangeable)
 expect_false(rule(list(sigma2 = "cauchy(0, 2)"), "gaussian")$exchangeable)
 expect_false(rule(list(phi1 = "gamma(2, 1)"), "gamma")$exchangeable)
 expect_false(rule(list(shape2 = "gamma(3, 1)"), "weibull", "survival")$exchangeable)
 expect_false(rule(list(scale1 = "gamma(3, 1)"), "weibull", "survival")$exchangeable)
 expect_true(rule(list(beta1 = "normal(0, 2.5)", beta2 = "normal(0,2.5)"), "binomial")$exchangeable)

 # Path A: a beta(a, b) prior with a != b identifies the labels (the starting
 # state is oriented towards the more probable labelling)
 r <- rule(NULL, "gaussian", m.rate = 0.6)
 expect_equal(r$orientation, "match-rate prior")
 expect_equal(r$pre_orient, "mass")
 expect_false(r$exchangeable)
 expect_equal(rule(NULL, "gaussian", m.rate = 0.2)$orientation, "match-rate prior")
 expect_equal(rule(list(theta = "beta(3, 7)"), "gaussian")$orientation, "match-rate prior")
 expect_true(rule(list(theta = "beta(2, 2)"), "gaussian")$exchangeable)
 expect_true(rule(NULL, "gaussian", m.rate = 0.5)$exchangeable)

 # Path B: a gamma_intercept or gamma_slope prior (or its alias gamma) with a
 # non-zero mean identifies the labels
 ruleB <- function(priors, ...) rule(priors, "gaussian", use_logistic = TRUE, M = 2L, ...)
 expect_true(ruleB(NULL)$exchangeable)
 expect_true(ruleB(list(gamma_intercept = "normal(0, 3)", gamma_slope = "normal(0, 1)"))$exchangeable)
 expect_equal(ruleB(list(gamma_slope = "normal(2, 1)"))$orientation, "match-rate prior")
 expect_equal(ruleB(list(gamma = "normal(-1, 1)"))$orientation, "match-rate prior")
 expect_equal(ruleB(list(gamma_intercept = "normal(0.5, 1)"))$orientation, "match-rate prior")
 expect_equal(ruleB(NULL, m.rate = 0.3)$orientation, "match-rate prior")

 # safe matches identify the labels whatever the priors
 r <- postlink:::.label_rule(engine_priors(NULL, "gaussian"), c(1L, integer(9)))
 expect_equal(r$orientation, "safe matches")
 expect_equal(r$pre_orient, "none")
 expect_false(r$exchangeable)
 r <- postlink:::.label_rule(engine_priors(NULL, "gaussian", m.rate = 0.6), c(1L, integer(9)))
 expect_equal(r$orientation, "safe matches")
})

test_that(".default_component_priors() tells the default component priors from the user's", {
 dflt <- function(priors, family, model_type = "glm", ...) {
  pf <- postlink:::prepare_mixbayes_priors(priors, family, model_type, ...)
  postlink:::.default_component_priors(pf, family, model_type)
 }
 expect_true(dflt(NULL, "binomial"))
 # a prior equal to the default, a theta prior, m.rate or linkage priors leave the
 # component priors at their defaults
 expect_true(dflt(list(beta1 = "normal(0, 2.5)"), "binomial"))
 expect_true(dflt(list(theta = "beta(2, 3)"), "binomial"))
 expect_true(dflt(NULL, "binomial", m.rate = 0.3))
 expect_true(dflt(list(gamma_slope = "normal(1, 1)"), "binomial", use_logistic = TRUE))
 expect_false(dflt(list(beta2 = "normal(0, 1)"), "binomial"))
 expect_false(dflt(list(intercept1 = "normal(0, 5)"), "binomial"))
 expect_true(dflt(NULL, "gaussian"))
 expect_false(dflt(list(sigma2 = "cauchy(0, 1)"), "gaussian"))
 expect_true(dflt(NULL, "gamma"))
 expect_false(dflt(list(phi1 = "gamma(2, 1)"), "gamma"))
 expect_true(dflt(NULL, "weibull", "survival"))
 expect_false(dflt(list(scale2 = "gamma(3, 1)"), "weibull", "survival"))
 expect_false(dflt(list(shape1 = "gamma(3, 1)"), "weibull", "survival"))
 expect_false(dflt(list(phi2 = "exponential(2)"), "gamma", "survival"))
})

# ------------------------------------------------------------------------------
# (a) An m.rate prior identifies the labels: no relabelling, so the posterior
# of theta is not cut at 0.5 (the ECR step used to exchange individual draws);
# the starting state is oriented towards the more probable labelling
# ------------------------------------------------------------------------------

test_that("with m.rate the draws are reported as sampled and theta is not truncated at 0.5", {
 LD <- as.data.frame(LD1000)
 set.seed(1)
 idx <- sample(nrow(LD), 400)
 X <- cbind("(Intercept)" = 1, BMI = as.numeric(scale(LD$BMI)), Age = as.numeric(scale(LD$Age)),
            Treatment = LD$Treatment)[idx, ]
 y <- LD$Disease_Status[idx]
 ctl <- list(iterations = 1500, burnin.iterations = 500, seed = 1)
 f <- suppressWarnings(glmMixBayes(X, y, "binomial", m.rate = 0.6, control = ctl))
 d <- f$diagnostics
 expect_equal(d$orientation, "match-rate prior")
 expect_false(d$relabelled)
 expect_identical(d$pre_oriented, isTRUE(d$orient_lp_change > 0))
 expect_false(d$orient_refused)
 expect_equal(d$n_swapped, 0L)
 expect_gte(d$ecr_share, 0)
 expect_true(d$other_labelling > 0 && d$other_labelling < 1)
 # the reported draws are the sampled draws
 pf <- postlink:::prepare_mixbayes_priors(NULL, "binomial", "glm", m.rate = 0.6)
 ep <- postlink:::build_engine_priors(pf, "binomial", "glm", K = ncol(X))
 post <- postlink:::run_mixbayes_engine("binomial", X, y, NULL, NULL, integer(nrow(X)), ep,
                                        postlink:::.mixbayes_settings(list(), ctl), pre_orient = "mass")
 expect_equal(f$estimates$theta, post$theta)
 expect_equal(unname(f$estimates$coefficients), unname(post$beta1))
 expect_identical(unname(f$m_samples), unname(post$z))
 expect_equal(d$lp, post$lp)
 # the prior puts 16% of its mass above 0.5 and the data barely move it: the
 # posterior is not cut at 0.5
 expect_gt(mean(f$estimates$theta > 0.5), 0.05)
 expect_gt(stats::quantile(f$estimates$theta, 0.975), 0.55)
})

# ------------------------------------------------------------------------------
# (b) Binomial defaults (component priors differ): oriented before sampling and
# relabelled after it; diagnostics$lp is the log posterior of the reported
# (relabelled) draws under fit$priors
# ------------------------------------------------------------------------------

test_that("binomial defaults: the starting state is oriented and the draws are relabelled", {
 set.seed(31)
 n <- 150
 x <- stats::rnorm(n)
 m <- stats::rbinom(n, 1, 0.6) == 1
 y <- stats::rbinom(n, 1, stats::plogis(ifelse(m, -1 + 2.5 * x, 1)))
 X <- cbind("(Intercept)" = 1, x = x)
 ctl <- list(iterations = 400, burnin.iterations = 200, seed = 1)

 # engine: the labels of the starting state are exchanged exactly when fewer
 # than half of the records are allocated to component 1 in it (on average
 # over the last sweeps of the selected pilot chain); the binomial
 # default priors differ only mildly between the components, so the estimated
 # log ratio of the posterior probabilities of the two labellings stays well
 # within the tolerance and the exchange is never refused
 ep <- engine_priors(NULL, "binomial")
 oriented <- logical(0)
 for (seed in 1:4) {
  mc <- postlink:::.mixbayes_settings(list(), list(iterations = 60, burnin.iterations = 30, seed = seed))
  post <- postlink:::run_mixbayes_engine("binomial", X, y, NULL, NULL, integer(n), ep, mc, pre_orient = TRUE)
  expect_identical(post$pre_oriented, post$start_share1 < 0.5)
  expect_false(post$orient_refused)
  if (post$pre_oriented) {
   expect_gt(post$orient_lp_change, -postlink:::.ORIENT_TOL)
  } else {
   expect_true(is.na(post$orient_lp_change))
  }
  oriented <- c(oriented, post$pre_oriented)
  # without the request the starting state is never oriented
  post0 <- postlink:::run_mixbayes_engine("binomial", X, y, NULL, NULL, integer(n), ep, mc)
  expect_false(post0$pre_oriented)
  expect_false(post0$orient_refused)
  expect_equal(post0$start_share1, post$start_share1)
 }
 expect_true(any(oriented))

 q <- with_warnings(glmMixBayes(X, y, "binomial", control = ctl))
 f <- q$value
 d <- f$diagnostics
 expect_equal(d$orientation, "majority component")
 expect_identical(d$pre_oriented, isTRUE(d$orient_lp_change > -postlink:::.ORIENT_TOL))
 expect_false(d$orient_refused)
 expect_length(d$pilot_lp, 5L)
 expect_true(d$relabelled)
 # the reported draws are aligned (ECR finds nothing left to exchange), and
 # component 1 is the majority component
 expect_equal(mean(postlink:::ecr_iterative_two(f$m_samples)), 0)
 expect_gt(mean(f$m_samples == 1L), 0.5)
 expect_false(any(grepl("meant to be the majority|visited both labellings", q$warnings)))
 # the stored priors are those of the components as sampled ...
 expect_equal(f$priors[["beta1"]], "normal(0, 2.5)")
 expect_equal(f$priors[["beta2"]], "normal(0, 5)")
 # ... and diagnostics$lp is the log posterior of the reported draws under
 # them: the one recorded by the sampler, plus its change on exchanging the
 # labels for the relabelled draws
 expect_equal(lp_reported(f, X, y), d$lp, tolerance = 1e-8)
 post <- postlink:::run_mixbayes_engine("binomial", X, y, NULL, NULL, integer(n), ep,
                                        postlink:::.mixbayes_settings(list(), ctl), pre_orient = "majority")
 swapped <- rowSums(unname(f$m_samples) != post$z) == n
 expect_equal(sum(swapped), d$n_swapped)
 expect_equal(d$lp, post$lp + ifelse(swapped, post$exchange_lp, 0))
 expect_true(all(post$exchange_lp != 0))
})

test_that("binomial defaults: a chain that switches labellings is relabelled (regression)", {
 # 70% correct matches with coefficients (-0.5, 2, -1) and mismatches with
 # probability 0.5: under the binomial defaults (beta1 normal(0, 2.5), beta2
 # normal(0, 5)) the chain of seed 4 spends about 40% of its draws in the
 # other labelling. The draws used to be reported as sampled (coefficient
 # means of about (1.1, 0.5, -0.1), theta 0.56, with warnings that they should
 # not be interpreted); relabelled, they recover the correct-match model
 set.seed(5)
 N <- 500
 x1 <- stats::rnorm(N); x2 <- stats::rnorm(N)
 zt <- ifelse(stats::runif(N) < 0.7, 1, 2)
 y <- ifelse(zt == 1, stats::rbinom(N, 1, stats::plogis(-0.5 + 2 * x1 - 1 * x2)), stats::rbinom(N, 1, 0.5))
 X <- cbind("(Intercept)" = 1, x1 = x1, x2 = x2)
 q <- with_warnings(glmMixBayes(X, y, "binomial", control = list(iterations = 3000, burnin.iterations = 1000,
                                                                 seed = 4)))
 f <- q$value
 d <- f$diagnostics
 expect_true(d$relabelled)
 expect_gt(d$ecr_share, 0.2)
 expect_gt(d$n_swapped, 0.2 * nrow(f$m_samples))
 expect_false(any(grepl("visited both labellings|meant to be the majority|Split R-hat", q$warnings)))
 b <- colMeans(f$estimates$coefficients)
 expect_lt(abs(b[["(Intercept)"]] + 0.5), 0.3)
 expect_lt(abs(b[["x1"]] - 2.2), 0.4)
 expect_lt(abs(b[["x2"]] + 1.2), 0.3)
 expect_gt(mean(f$estimates$theta), 0.75)
 expect_gt(mean((colMeans(f$m_samples == 1L) > 0.5) == (zt == 1)), 0.7)
 # the log posterior of the reported draws is consistent with the relabelling
 expect_equal(lp_reported(f, X, y), d$lp, tolerance = 1e-8)
})

test_that("starting values are not reoriented and safe matches are never reoriented", {
 set.seed(9)
 n <- 200
 x <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, 0.3) == 1
 y <- ifelse(match, 1 + 2 * x + stats::rnorm(n, sd = 0.5), -3 - x + stats::rnorm(n, sd = 0.5))
 X <- cbind("(Intercept)" = 1, x = x)
 ep <- engine_priors(NULL, "gaussian")
 # starting values with component 1 on the 30% minority
 init <- list(beta1 = c(1, 2), beta2 = c(-3, -1), theta = 0.3, disp1 = 0.5, disp2 = 0.5)
 mc <- postlink:::.mixbayes_settings(list(), list(iterations = 60, burnin.iterations = 30, seed = 2, init = init))
 post <- postlink:::run_mixbayes_engine("gaussian", X, y, NULL, NULL, integer(n), ep, mc, pre_orient = TRUE)
 expect_lt(post$start_share1, 0.5)
 expect_false(post$pre_oriented)
 expect_false(post$orient_refused)
 expect_true(is.na(post$orient_lp_change))
 # with component priors that differ the draws are relabelled after sampling
 # all the same (majority convention): the global flip exchanges the labels of
 # every draw, and diagnostics$lp is the log posterior of the relabelled draws
 # (the priors of the two components are exchanged with them)
 q <- with_warnings(glmMixBayes(X, y, "gaussian", priors = list(beta2 = "normal(0, 4)"),
                                control = list(iterations = 400, burnin.iterations = 200, seed = 5, init = init)))
 d <- q$value$diagnostics
 expect_false(d$pre_oriented)
 expect_true(d$relabelled)
 expect_equal(d$n_swapped, nrow(q$value$m_samples))
 expect_gt(mean(q$value$m_samples == 1L), 0.5)
 expect_false(any(grepl("meant to be the majority", q$warnings)))
 expect_equal(lp_reported(q$value, X, y), d$lp, tolerance = 1e-8)
 # safe matches are fixed in component 1: never reoriented
 safe <- as.integer(match & seq_len(n) %% 4 == 0)
 mc <- postlink:::.mixbayes_settings(list(), list(iterations = 60, burnin.iterations = 30, seed = 2))
 post <- postlink:::run_mixbayes_engine("gaussian", X, y, NULL, NULL, safe, ep, mc, pre_orient = TRUE)
 expect_lt(post$start_share1, 0.5)
 expect_false(post$pre_oriented)
 expect_true(all(post$z[, safe == 1L] == 1L))
})

# ------------------------------------------------------------------------------
# (c) Gaussian defaults (exchangeable): ECR alignment and majority flip after
# sampling, which leave the log posterior of every draw unchanged
# ------------------------------------------------------------------------------

test_that("exchangeable priors: draws are relabelled after sampling and lp is invariant", {
 set.seed(9)
 n <- 200
 x <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, 0.3) == 1
 y <- ifelse(match, 1 + 2 * x + stats::rnorm(n, sd = 0.5), -3 - x + stats::rnorm(n, sd = 0.5))
 X <- cbind("(Intercept)" = 1, x = x)
 ctl <- list(iterations = 400, burnin.iterations = 200, seed = 5)
 # starting values put component 1 on the 30% minority: they are respected
 # (no orientation before sampling), and the global flip after sampling
 # exchanges the labels of every draw
 init <- list(beta1 = c(1, 2), beta2 = c(-3, -1), theta = 0.3, disp1 = 0.5, disp2 = 0.5)
 utils::capture.output(
  expect_message(f <- glmMixBayes(X, y, "gaussian", control = c(ctl, list(init = init, verbose = TRUE))),
                 "Global label swap"))
 d <- f$diagnostics
 expect_equal(d$orientation, "majority component")
 expect_false(d$pre_oriented)
 expect_true(d$relabelled)
 expect_equal(d$n_swapped, nrow(f$m_samples))
 expect_equal(d$ecr_share, 0)
 expect_gt(mean(f$estimates$theta), 0.6)
 expect_lt(abs(mean(f$estimates$coefficients[, "x"]) + 1), 0.2)    # the majority line -3 - x
 expect_equal(f$priors[["beta1"]], f$priors[["beta2"]])
 # the log posterior recorded for the sampled draws is that of the relabelled draws
 expect_equal(lp_reported(f, X, y), d$lp, tolerance = 1e-8)

 # without starting values: the same rules, and the same invariance
 f2 <- suppressMessages(glmMixBayes(X, y, "gaussian", control = ctl))
 expect_true(f2$diagnostics$relabelled)
 expect_gt(mean(f2$m_samples == 1L), 0.5)
 expect_equal(lp_reported(f2, X, y), f2$diagnostics$lp, tolerance = 1e-8)
})

# ------------------------------------------------------------------------------
# (d) A gamma_slope prior with a non-zero mean identifies the labels: no
# majority flip (it used to report the mismatch component as component 1)
# ------------------------------------------------------------------------------

test_that("a gamma_slope prior with a non-zero mean identifies the components", {
 set.seed(101)
 N <- 300
 match <- stats::runif(N) < 0.4                   # correct matches are the minority
 s <- ifelse(match, stats::rnorm(N, 1, 1), stats::rnorm(N, -1, 1))
 x <- stats::rnorm(N)
 y <- 1 + 2 * x + stats::rnorm(N, sd = 0.5)
 y[!match] <- sample(y[!match])
 X <- cbind("(Intercept)" = 1, x = x)
 Z <- cbind("(Intercept)" = 1, s = s)
 ctl <- list(iterations = 1000, burnin.iterations = 500, seed = 1)
 for (pr in list(list(gamma_slope = "normal(2, 1)"), list(gamma = "normal(2, 1)"))) {
  f <- suppressWarnings(glmMixBayes(X, y, "gaussian", Z = Z, priors = pr, control = ctl))
  d <- f$diagnostics
  expect_equal(d$orientation, "match-rate prior")
  expect_identical(d$pre_oriented, isTRUE(d$orient_lp_change > 0))
  expect_false(d$relabelled)
  expect_equal(d$n_swapped, 0L)
  # the slope prior identifies the labels firmly here
  expect_lt(d$other_labelling, 0.05)
  # component 1 holds the (minority) correct matches, with a positive slope
  # on the linkage score as the prior says
  expect_gt(mean(f$estimates$gamma[, "s"]), 1)
  expect_gt(mean(f$m_samples[, match] == 1L), 0.7)
  expect_lt(mean(f$m_samples[, !match] == 1L), 0.3)
  expect_lt(mean(f$m_samples == 1L), 0.5)
 }
 # without the slope prior the majority convention applies
 f0 <- suppressWarnings(glmMixBayes(X, y, "gaussian", Z = Z, control = ctl))
 expect_equal(f0$diagnostics$orientation, "majority component")
 expect_true(f0$diagnostics$relabelled)
 expect_gt(mean(f0$m_samples == 1L), 0.5)
})

# ------------------------------------------------------------------------------
# (e) Safe matches: no relabelling; the ECR diagnostic warns only above 10%
# ------------------------------------------------------------------------------

test_that("safe matches: no relabelling, and a warning only when the ECR share exceeds 0.1", {
 ctl <- list(iterations = 3000, burnin.iterations = 1000, seed = 1)
 # two groups with opposite slopes and a single safe match: a genuinely
 # bimodal posterior over the other records
 set.seed(1)
 n <- 30
 x <- stats::rnorm(n)
 grp <- rep(c("A", "B"), each = n / 2)
 y <- ifelse(grp == "A", 0.5 + x, -0.5 - x) + stats::rnorm(n)
 safe <- integer(n)
 safe[which.min(abs(x[grp == "A"]))] <- 1L
 q <- with_warnings(glmMixBayes(cbind(1, x), y, "gaussian", safe.matches = safe, control = ctl))
 d <- q$value$diagnostics
 expect_equal(d$orientation, "safe matches")
 expect_false(d$relabelled)
 expect_false(d$pre_oriented)
 expect_equal(d$n_swapped, 0L)
 expect_equal(d$ecr_share, mean(postlink:::ecr_iterative_two(q$value$m_samples, cols = which(safe == 0L))))
 expect_gt(d$ecr_share, 0.1)
 expect_true(any(grepl(sprintf("visited both labellings of the mixture: %.1f%%", 100 * d$ecr_share),
                       q$warnings, fixed = TRUE)))
 expect_true(all(q$value$m_samples[, safe == 1L] == 1L))
 # summary() notes that the reported draws mix both labellings
 sm <- summary(q$value)
 expect_equal(sm$mixed_labels, d$ecr_share)
 expect_output(print(sm), sprintf("visited both labellings of the mixture (%.1f%%", 100 * d$ecr_share),
               fixed = TRUE)

 # well separated components: no draw is exchanged and nothing is reported
 set.seed(2)
 n <- 120
 x <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, 0.75) == 1
 y <- ifelse(match, 1 + 2 * x, -3 - x) + stats::rnorm(n, sd = 0.4)
 safe <- as.integer(match & seq_len(n) %% 5 == 0)
 q <- with_warnings(glmMixBayes(cbind(1, x), y, "gaussian", safe.matches = safe,
                                control = list(iterations = 600, burnin.iterations = 300, seed = 3)))
 expect_lte(q$value$diagnostics$ecr_share, 0.1)
 expect_false(any(grepl("visited both labellings", q$warnings)))
 expect_null(summary(q$value)$mixed_labels)
 expect_false(any(grepl("visited both labellings", utils::capture.output(print(summary(q$value))))))
 expect_equal(q$value$diagnostics$n_swapped, 0L)
})

test_that("the ECR share and its warning use the records not flagged as safe matches", {
 S <- 20; N <- 12
 safe <- c(rep(1L, 6), rep(0L, 6))
 z <- matrix(c(rep(1L, 6), 1L, 1L, 1L, 2L, 2L, 2L), S, N, byrow = TRUE)
 z[1:4, 7:12] <- 3L - z[1:4, 7:12]                    # 4 of 20 draws: 20%
 expect_equal(mean(postlink:::ecr_iterative_two(z, cols = 7:12)), 0.2)
 expect_equal(postlink:::ecr_iterative_two(z, cols = 7:12),
              postlink:::ecr_iterative_two(z[, 7:12]))
 expect_error(postlink:::ecr_iterative_two(z, cols = 13L), "out of range")
 expect_warning(out <- postlink:::align_mixture_labels(z, list(), theta = rep(0.5, S), safe = safe),
                "20.0%")
 expect_equal(out$ecr_share, 0.2)
 # the warning names the help page of the fitting function
 expect_warning(postlink:::align_mixture_labels(z, list(), safe = safe), "in ?glmMixBayes)", fixed = TRUE)
 expect_warning(postlink:::align_mixture_labels(z, list(), safe = safe, help = "survregMixBayes"),
                "in ?survregMixBayes)", fixed = TRUE)
 expect_identical(out$z, z)
 # the same draws with an informative prior on the match probability, and
 # with component priors that differ (majority convention, no relabelling)
 expect_warning(postlink:::align_mixture_labels(z[, 7:12], list(), orientation = "match-rate prior",
                                                relabel = FALSE),
                "prior on the match probability identifies component 1")
 expect_warning(out2 <- postlink:::align_mixture_labels(z[, 7:12], list(), orientation = "majority component",
                                                        relabel = FALSE),
                "component-specific priors identify component 1")
 expect_false(out2$relabelled)
 expect_false(any(out2$swapped))
 # a share of at most 10% is not reported
 z1 <- z; z1[3:4, 7:12] <- 3L - z1[3:4, 7:12]         # 2 of 20 draws: 10%
 expect_silent(out3 <- postlink:::align_mixture_labels(z1, list(), safe = safe))
 expect_equal(out3$ecr_share, 0.1)
})

test_that("survregMixBayes() applies the same rules and records them", {
 set.seed(4)
 n <- 160
 x <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, 0.75) == 1
 tt <- stats::rweibull(n, shape = 1.8, scale = exp(ifelse(match, 0.5 + 0.8 * x, 2)))
 X <- cbind("(Intercept)" = 1, x = x)
 ctl <- list(iterations = 400, burnin.iterations = 200, seed = 2)
 f <- suppressWarnings(survregMixBayes(X, cbind(tt, 1L), dist = "weibull", control = ctl))
 d <- f$diagnostics
 expect_equal(d$orientation, "majority component")
 expect_true(d$relabelled)
 expect_true(is.logical(d$pre_oriented))
 expect_true(d$ecr_share >= 0 && d$ecr_share <= 1)
 expect_gt(mean(f$m_samples == 1L), 0.5)
 fm <- suppressWarnings(survregMixBayes(X, cbind(tt, 1L), dist = "weibull", m.rate = 0.3, control = ctl))
 expect_equal(fm$diagnostics$orientation, "match-rate prior")
 expect_false(fm$diagnostics$relabelled)
 expect_equal(fm$diagnostics$n_swapped, 0L)
 fs <- suppressWarnings(survregMixBayes(X, cbind(tt, 1L), dist = "weibull", priors = list(scale2 = "gamma(3, 1)"),
                                        control = ctl))
 expect_equal(fs$diagnostics$orientation, "majority component")
 # component priors that differ: relabelled after sampling unless the
 # exchange of the starting state was refused
 expect_false(fs$diagnostics$orient_refused)
 expect_true(fs$diagnostics$relabelled)
 expect_gt(mean(fs$m_samples == 1L), 0.5)
 expect_length(fs$diagnostics$pilot_lp, 5L)
 # summary() notes draws that mix both labellings (only when not relabelled)
 fs$diagnostics$ecr_share <- 0.25
 expect_null(summary(fs)$mixed_labels)
 fs$diagnostics$relabelled <- FALSE
 expect_output(print(summary(fs)), "visited both labellings of the mixture (25.0%", fixed = TRUE)
 expect_output(print(summary(fs)), "see 'Label switching' in ?survregMixBayes).", fixed = TRUE)
 # the note names the help page it is given and ends with a newline
 out <- utils::capture.output({
  postlink:::.print_mixed_labels(0.25, "survregMixBayes")
  cat("next\n")
 })
 expect_length(out, 2L)
 expect_match(out[1], "see 'Label switching' in ?survregMixBayes).", fixed = TRUE)
 expect_identical(out[2], "next")
 expect_output(postlink:::.print_mixed_labels(0.25), "?glmMixBayes)", fixed = TRUE)
 expect_silent(postlink:::.print_mixed_labels(NULL))
 expect_null(postlink:::.mixed_labels(list(relabelled = TRUE, ecr_share = 0.25)))
 expect_null(postlink:::.mixed_labels(list(relabelled = FALSE, ecr_share = 0.1)))
 expect_null(postlink:::.mixed_labels(list()))

 # a tight beta2 prior (mismatches without slope) and a 30% minority of
 # correct matches with a slope: the exchange of the starting state is refused
 # (the check chains confirm it) and the warnings name ?survregMixBayes
 set.seed(4)
 N <- 200
 xs <- stats::rnorm(N)
 ms <- stats::runif(N) < 0.3
 tt2 <- ifelse(ms, stats::rweibull(N, 2, exp(1 + 1.5 * xs)), stats::rweibull(N, 2, exp(1)))
 q <- with_warnings(survregMixBayes(cbind("(Intercept)" = 1, x = xs), cbind(tt2, 1L), dist = "weibull",
                                    priors = list(beta2 = "normal(0, 0.01)"),
                                    control = list(seed = 1, iterations = 400, burnin.iterations = 200)))
 d <- q$value$diagnostics
 expect_true(d$orient_refused)
 expect_false(d$pre_oriented)
 expect_true(is.finite(d$orient_check[["kept"]]))
 expect_gt(mean((colMeans(q$value$m_samples == 1L) > 0.5) == ms), 0.75)
 expect_true(any(grepl("favour the labelling in which component 1 is the minority", q$warnings)))
 expect_true(any(grepl("in ?survregMixBayes)", q$warnings, fixed = TRUE)))
 expect_false(any(grepl("?glmMixBayes", q$warnings, fixed = TRUE)))
})

# ------------------------------------------------------------------------------
# (f) Orienting the starting state keeps the sampler valid
# ------------------------------------------------------------------------------

test_that("the orientation of the starting state leaves the posterior unchanged", {
 # component 2 is pinned at N(5, 1) and component 1 holds about 30% of the
 # records, so every pilot state has component 2 as the majority and its labels
 # are exchanged before the main chain (orient_tol = Inf forces the exchange,
 # which the default tolerance refuses under these pinned priors); the chain
 # must return to the posterior, which (checked against the exact posterior by
 # quadrature in the review of the sampler) is the one sampled without the
 # orientation
 set.seed(202)
 N <- 120
 x <- stats::rnorm(N)
 zt <- ifelse(stats::runif(N) < 0.3, 1, 2)
 y <- ifelse(zt == 1, 1 + 0.5 * x + stats::rnorm(N, 0, 0.7), 5 + stats::rnorm(N, 0, 1))
 X <- cbind(1, x)
 pin <- 1e-3
 ep <- list(beta1_mean = c(0, 0), beta1_sd = c(10, 10), beta2_mean = c(5, 0), beta2_sd = c(pin, pin),
            theta = c(2, 1.5), disp1 = c(2, 0, 2.5, 0), disp2 = c(5, 0, pin, 0), intercept = TRUE)
 draws <- function(pre_orient, seed, tol = Inf) {
  mc <- postlink:::.mixbayes_settings(list(), list(iterations = 12000, burnin.iterations = 2000, seed = seed))
  p <- postlink:::run_mixbayes_engine("gaussian", X, y, NULL, NULL, integer(N), ep, mc, pre_orient = pre_orient,
                                      orient_tol = tol)
  list(pre = p$pre_oriented, refused = p$orient_refused, check = p$orient_check,
       d = cbind(b0 = p$beta1[, 1], b1 = p$beta1[, 2], sigma1 = p$disp1, theta = p$theta))
 }
 a <- draws(TRUE, 11)
 b <- draws(FALSE, 12)
 expect_true(a$pre)
 expect_false(b$pre)
 # with orient_tol = Inf (and without orientation) no check chains are run
 expect_true(all(is.na(a$check)) && all(is.na(b$check)))
 mcse <- function(v) stats::sd(v) / sqrt(postlink:::.ess_ar(v))
 for (j in colnames(a$d)) {
  tol <- 4 * sqrt(mcse(a$d[, j])^2 + mcse(b$d[, j])^2)
  expect_lt(abs(mean(a$d[, j]) - mean(b$d[, j])), tol, label = j)
  expect_lt(abs(stats::sd(a$d[, j]) / stats::sd(b$d[, j]) - 1), 0.1, label = j)
 }
 # exact posterior means (quadrature; review of the sampler, 4 x 60000 draws agree)
 expect_lt(abs(mean(a$d[, "theta"]) - 0.3203), 4 * mcse(a$d[, "theta"]))
 expect_lt(abs(mean(a$d[, "b0"]) - 1.0709), 4 * mcse(a$d[, "b0"]))
 # with the default tolerance the estimated change is far below -2, so the two
 # check chains are run before the main chain (they consume random numbers and
 # leave the sampler state behind); they refuse the exchange, and the main
 # chain still samples the exact posterior
 cc <- draws(TRUE, 13, tol = postlink:::.ORIENT_TOL)
 expect_true(cc$refused)
 expect_false(cc$pre)
 expect_named(cc$check, c("exchanged", "kept", "margin", "share1"))
 expect_true(all(is.finite(cc$check[c("kept", "margin", "share1")])))
 expect_true(cc$check[["share1"]] < 0.5 || !(cc$check[["exchanged"]] >= cc$check[["kept"]] - cc$check[["margin"]]))
 expect_lt(abs(mean(cc$d[, "theta"]) - 0.3203), 4 * mcse(cc$d[, "theta"]))
 expect_lt(abs(mean(cc$d[, "b0"]) - 1.0709), 4 * mcse(cc$d[, "b0"]))
})

# ------------------------------------------------------------------------------
# (g) The orientation is refused when the component-specific priors make the
# labelling with component 1 as the majority much less probable
# ------------------------------------------------------------------------------

test_that("the starting state is not oriented against strongly asymmetric component priors", {
 # 30% of the records follow y = 1 + 2x, 70% are flat; beta2 = normal(0, 0.01)
 # says that component 2 has no slope, so the flat majority can only be
 # component 2. Exchanging the labels of a pilot state used to start the main
 # chain in a mirror mode about 65 log-units lower that it did not leave
 # within 10,000 iterations (the slope-2 group was lost, silently)
 set.seed(9)
 n <- 200
 x <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, 0.3) == 1
 X <- cbind("(Intercept)" = 1, x = x)
 y0 <- ifelse(match, 1 + 2 * x + stats::rnorm(n, sd = 0.5), 0.5 + stats::rnorm(n, sd = 0.5))
 eng <- function(priors, seed, tol = postlink:::.ORIENT_TOL, it = 1500) {
  mc <- postlink:::.mixbayes_settings(list(), list(iterations = it, burnin.iterations = 500, seed = seed))
  postlink:::run_mixbayes_engine("gaussian", X, y0, NULL, NULL, integer(n), engine_priors(priors, "gaussian"),
                                 mc, pre_orient = TRUE, orient_tol = tol)
 }
 for (pr in list(list(beta2 = "normal(0, 0.01)"), list(beta2 = "normal(0, 0.1)"))) {
  for (seed in 1:3) {
   p <- eng(pr, seed)
   expect_lt(p$start_share1, 0.5)
   expect_true(p$orient_refused)
   expect_false(p$pre_oriented)
   expect_lt(p$orient_lp_change, -100)
   # the check chains confirm the refusal: the chain started with the labels
   # exchanged returned to the original labelling, or stayed far below the
   # chain started without exchange
   chk <- p$orient_check
   expect_true(chk[["margin"]] >= postlink:::.ORIENT_TOL)
   expect_true(chk[["share1"]] < 0.5 || chk[["exchanged"]] < chk[["kept"]] - chk[["margin"]])
   # the main chain stays at the level of the selected pilot chain and keeps
   # the slope-2 group in component 1
   expect_gt(mean(p$lp), max(p$pilot_lp) - 5)
   expect_lt(abs(mean(p$beta1[, 2]) - 2), 0.2)
   expect_gt(mean(p$z[, match] == 1L), 0.7)
   expect_lt(mean(p$z[, !match] == 1L), 0.25)
  }
 }
 # identical component priors: the estimated change is exactly 0, so the
 # exchange is always made (and undone or confirmed by relabelling after sampling)
 for (seed in 1:3) {
  p <- eng(NULL, seed, it = 600)
  if (p$start_share1 < 0.5) {
   expect_identical(p$orient_lp_change, 0)
   expect_true(p$pre_oriented)
  } else {
   expect_true(is.na(p$orient_lp_change))
  }
  expect_false(p$orient_refused)
 }
 # weakly asymmetric priors: oriented at once, without check chains
 p <- eng(list(beta1 = "normal(0, 3)", beta2 = "normal(0, 4)"), 1, it = 600)
 expect_lt(p$start_share1, 0.5)
 expect_true(p$pre_oriented)
 expect_gt(p$orient_lp_change, -postlink:::.ORIENT_TOL)
 expect_true(all(is.na(p$orient_check)))
 # beta2 = normal(0, 0.5): moving the slope-2 group to component 2 lowers the
 # log posterior by about 8, which the estimate and the check chains (which mix
 # well here, so that the margin is the tolerance itself) agree on
 for (seed in 1:2) {
  p <- eng(list(beta2 = "normal(0, 0.5)"), seed, it = 600)
  expect_true(p$orient_refused)
  expect_lt(p$orient_lp_change, -5)
  chk <- p$orient_check
  expect_gt(chk[["share1"]], 0.5)
  expect_lt(chk[["exchanged"]] - chk[["kept"]], -5)
  expect_lt(chk[["margin"]], 4)
 }

 # fit level: the draws are reported as sampled, the slope-2 group is
 # component 1, and the warning on the minority component explains why
 q <- with_warnings(glmMixBayes(X, y0, "gaussian", priors = list(beta2 = "normal(0, 0.01)"),
                                control = list(iterations = 1500, burnin.iterations = 500, seed = 2)))
 f <- q$value
 d <- f$diagnostics
 expect_equal(d$orientation, "majority component")
 expect_true(d$orient_refused)
 expect_false(d$pre_oriented)
 # the draws are aligned, but not exchanged globally either: the reported
 # labelling holds practically all of the posterior probability
 expect_true(d$relabelled)
 expect_equal(d$n_swapped, 0L)
 expect_lt(d$other_labelling, 0.01)
 expect_lt(d$orient_lp_change, -100)
 expect_gt(mean(d$lp), max(d$pilot_lp) - 5)
 expect_lt(abs(mean(f$estimates$coefficients[, "x"]) - 2), 0.2)
 expect_gt(mean((colMeans(f$m_samples == 1L) > 0.5) == match), 0.8)
 expect_true(any(grepl("favour the labelling in which component 1 is the minority", q$warnings)))
 expect_false(any(grepl("started with component 1 as the majority", q$warnings)))
 # with the expected mismatch rate the match-rate prior identifies the
 # components: the starting state is oriented towards the more probable
 # labelling (never refused), and no warning on the minority component or on
 # a weak identification is issued
 q2 <- with_warnings(glmMixBayes(X, y0, "gaussian", priors = list(beta2 = "normal(0, 0.01)"), m.rate = 0.7,
                                 control = list(iterations = 1500, burnin.iterations = 500, seed = 2)))
 expect_equal(q2$value$diagnostics$orientation, "match-rate prior")
 expect_false(q2$value$diagnostics$orient_refused)
 expect_identical(q2$value$diagnostics$pre_oriented, q2$value$diagnostics$orient_lp_change > 0)
 expect_lt(q2$value$diagnostics$other_labelling, 0.05)
 expect_false(any(grepl("majority convention|only weakly|less probable", q2$warnings)))
 expect_lt(abs(mean(q2$value$estimates$coefficients[, "x"]) - 2), 0.2)
})

test_that("binomial defaults with many covariates: component 1 is the majority for every seed", {
 # 15 covariates with slopes of +-1 for the 70% correct matches. Under the
 # binomial default priors (beta1 normal(0, 2.5), beta2 normal(0, 5)) the
 # estimated log ratio of the two labellings falls far below -2 for selected
 # pilot states with component 2 as the majority (its slopes are inflated to
 # about 4), but a chain started with the labels exchanged moves on to a mode
 # of comparable posterior probability: the check chains accept the exchange,
 # so the reported correct-match component does not depend on the seed
 set.seed(1503)
 n <- 400
 K <- 15
 Xc <- matrix(stats::rnorm(n * K), n, K, dimnames = list(NULL, paste0("x", seq_len(K))))
 m <- stats::runif(n) < 0.7
 b <- rep(c(1, -1), length.out = K)
 y <- stats::rbinom(n, 1, stats::plogis(ifelse(m, drop(Xc %*% b), 0.3)))
 X <- cbind("(Intercept)" = 1, Xc)
 checked <- logical(0)
 for (seed in 3:4) {
  q <- with_warnings(glmMixBayes(X, y, "binomial", control = list(iterations = 600, burnin.iterations = 300,
                                                                  seed = seed)))
  f <- q$value
  d <- f$diagnostics
  expect_equal(d$orientation, "majority component")
  expect_true(d$relabelled)
  expect_false(d$orient_refused)
  expect_gt(mean(f$m_samples == 1L), 0.6)
  expect_gt(mean((colMeans(f$m_samples == 1L) > 0.5) == m), 0.65)
  # component 1 carries the slopes of the correct matches
  expect_gt(mean(f$estimates$coefficients[, "x1"]), 0.5)
  expect_lt(mean(f$estimates$coefficients[, "x2"]), -0.5)
  expect_false(any(grepl("meant to be the majority", q$warnings)))
  expect_false(any(grepl("trapped in a minor mode", q$warnings)))
  if (isTRUE(d$pre_oriented) && isTRUE(d$orient_lp_change < -postlink:::.ORIENT_TOL)) {
   chk <- d$orient_check
   expect_gt(chk[["share1"]], 0.5)
   expect_gte(chk[["exchanged"]], chk[["kept"]] - chk[["margin"]])
   checked <- c(checked, TRUE)
  }
 }
 expect_true(any(checked))
})

test_that("the warning on a minority component 1 is worded after the starting state", {
 msg <- function(start) postlink:::.minority_message(0.3, start)
 expect_match(msg(list(refused = TRUE, lp_change = -8.04, share1 = 0.32)),
              "lower by about 8.0 (more than 2)", fixed = TRUE)
 expect_match(msg(list(init = TRUE, share1 = 0.3, pre_oriented = FALSE, refused = FALSE)),
              "starting values in control$init put component 1 on the minority", fixed = TRUE)
 # starting values with component 1 as the majority: the chain moved
 expect_match(msg(list(init = TRUE, share1 = 0.7, pre_oriented = FALSE, refused = FALSE)),
              "started with component 1 as the majority")
 expect_match(msg(list(init = FALSE, share1 = 0.3, pre_oriented = TRUE, refused = FALSE)),
              "started with component 1 as the majority")
 neutral <- msg(NULL)
 expect_match(neutral, "^Only 30% of the records")
 expect_false(grepl("started with|starting values|favour", neutral))
 # refusals: the result of the check chains, the help page, an infinite change,
 # and the default priors of the family
 refused <- list(refused = TRUE, lp_change = -8.04, share1 = 0.32,
                 check = c(exchanged = -231, kept = -223, margin = 2, share1 = 0.7))
 m1 <- msg(refused)
 expect_match(m1, "stayed 8.0 below one started without them exchanged (more than the margin of 2.0)",
              fixed = TRUE)
 expect_match(m1, "?glmMixBayes)", fixed = TRUE)
 expect_match(m1, "If these priors describe the two components as intended", fixed = TRUE)
 m2 <- postlink:::.minority_message(0.3, modifyList(refused, list(check = c(exchanged = -223, kept = -223,
                                                                             margin = 2, share1 = 0.3))),
                                    help = "survregMixBayes")
 expect_match(m2, "returned to the labelling in which component 1 is the minority", fixed = TRUE)
 expect_match(m2, "?survregMixBayes)", fixed = TRUE)
 expect_match(msg(modifyList(refused, list(check = c(exchanged = -Inf, kept = -223, margin = 2, share1 = 0.7)))),
              "did not reach a finite log posterior", fixed = TRUE)
 m3 <- msg(list(refused = TRUE, lp_change = -Inf, share1 = 0.32))
 expect_match(m3, "infinitely lower", fixed = TRUE)
 expect_false(grepl("Inf|check chain", m3))
 m4 <- msg(modifyList(refused, list(default_priors = TRUE)))
 expect_match(m4, "These are the default priors", fixed = TRUE)
 expect_match(m4, "beta1 = normal(0, 2.5) and beta2 = normal(0, 5)", fixed = TRUE)
 expect_match(m4, "'m.rate'", fixed = TRUE)
 expect_false(grepl("as intended", m4))
})

test_that("a main chain far below the selected pilot chain after the orientation is reported", {
 post <- list(pre_oriented = TRUE, pilot_lp = c(-101, -100, -130), lp = rep(-112, 20), orient_lp_change = -1)
 expect_warning(expect_true(postlink:::.check_orientation_lp(post, exchangeable = FALSE)),
                "trapped in a minor mode")
 # within the margin (the expected level includes a negative estimated change)
 post$lp <- rep(-105.5, 20)
 expect_silent(expect_false(postlink:::.check_orientation_lp(post, exchangeable = FALSE)))
 post$lp <- rep(-112, 20)
 # exchangeable components (the exchange leaves the posterior unchanged), an
 # orientation towards the more probable labelling under the match-rate
 # prior, not oriented, or no pilot chains: no check
 expect_silent(expect_false(postlink:::.check_orientation_lp(post, exchangeable = TRUE)))
 expect_silent(expect_false(postlink:::.check_orientation_lp(post, FALSE, orientation = "match-rate prior")))
 expect_silent(postlink:::.check_orientation_lp(modifyList(post, list(pre_oriented = FALSE)), FALSE))
 expect_silent(postlink:::.check_orientation_lp(modifyList(post, list(pilot_lp = numeric(0))), FALSE))
 # the margin grows with the standard deviation of the log posterior of the
 # stored draws: a gap of 9.5 is within twice an SD of 8 (large logistic
 # models), but not within the minimum margin of 5 for an SD of 1
 lp_sd <- function(m, s) m + s * (rep(c(-1, 1), 50))
 post <- list(pre_oriented = TRUE, pilot_lp = c(-1330, -1340), orient_lp_change = 11)
 expect_silent(expect_false(postlink:::.check_orientation_lp(modifyList(post, list(lp = lp_sd(-1339.5, 8))), FALSE)))
 expect_warning(expect_true(postlink:::.check_orientation_lp(modifyList(post, list(lp = lp_sd(-1339.5, 1))), FALSE)),
                "more than the margin of 5.0")
 # after check chains, the expected level allows for their difference
 post <- list(pre_oriented = TRUE, pilot_lp = -100, lp = rep(-108, 20), orient_lp_change = -12,
              orient_check = c(exchanged = -104, kept = -100, margin = 6, share1 = 0.8))
 expect_silent(expect_false(postlink:::.check_orientation_lp(post, FALSE)))
 post$orient_check[["exchanged"]] <- -100
 expect_warning(postlink:::.check_orientation_lp(post, FALSE), "trapped in a minor mode")
 # without (finite) check results, the estimated change is used
 post$orient_check <- c(exchanged = NA, kept = NA, margin = NA, share1 = NA)
 expect_silent(postlink:::.check_orientation_lp(post, FALSE))
})

# ------------------------------------------------------------------------------
# (h) Under the match-rate prior the starting state is oriented towards the
# more probable labelling, and a weak identification is reported
# ------------------------------------------------------------------------------

test_that("the log posterior of relabelled draws and the share of the other labelling", {
 lp <- c(-10, -11, -12)
 ex <- c(-2, 1, 0)
 r <- postlink:::.relabel_lp(lp, ex, c(FALSE, TRUE, FALSE))
 expect_equal(r$lp, c(-10, -10, -12))
 expect_equal(r$exchange_lp, c(-2, -1, 0))
 # nothing swapped, or no exchange changes (safe matches): unchanged
 expect_identical(postlink:::.relabel_lp(lp, ex, logical(3))$lp, lp)
 expect_identical(postlink:::.relabel_lp(lp, NULL, rep(TRUE, 3))$lp, lp)
 # log mean exp of the exchange changes, as a share of the two labellings
 expect_equal(postlink:::.other_labelling_share(c(0, 0)), 0.5)
 expect_equal(postlink:::.other_labelling_share(c(log(0.1), log(0.3))), 0.2 / 1.2)
 expect_true(is.na(postlink:::.other_labelling_share(rep(NA_real_, 3))))
 expect_true(is.na(postlink:::.other_labelling_share(NULL)))
 expect_equal(postlink:::.other_labelling_share(c(-Inf, -Inf)), 0)
 expect_equal(postlink:::.other_labelling_share(c(-800, -801)), 0)
 expect_equal(postlink:::.other_labelling_share(c(800, 700)), 1)

 # the warning under the match-rate prior
 expect_silent(expect_false(postlink:::.check_label_mass(0.04, "match-rate prior")))
 expect_warning(expect_true(postlink:::.check_label_mass(0.2, "match-rate prior")),
                "identifies the components only weakly: the labelling with the two components exchanged holds about 20%")
 expect_warning(postlink:::.check_label_mass(0.7, "match-rate prior", help = "survregMixBayes"),
                "less probable of the two labellings.*component 1 may describe.*[?]survregMixBayes")
 expect_warning(postlink:::.check_label_mass(0.95, "match-rate prior"), "component 1 probably describes the mismatches")
 # not under the other rules, without an estimate, or when the draws mix both
 # labellings (reported by the ECR warning)
 expect_silent(postlink:::.check_label_mass(0.5, "majority component"))
 expect_silent(postlink:::.check_label_mass(0.4, "safe matches"))
 expect_silent(postlink:::.check_label_mass(NA_real_, "match-rate prior"))
 expect_silent(postlink:::.check_label_mass(0.4, "match-rate prior", ecr_share = 0.2))

 # log mean exp
 expect_equal(postlink:::.log_mean_exp(c(log(2), log(4), NA)), log(3))
 expect_true(is.na(postlink:::.log_mean_exp(c(NA, NA))))
 expect_equal(postlink:::.log_mean_exp(c(-Inf, -Inf)), -Inf)
})

test_that("the global exchange after sampling is refused when the component priors identify the labels", {
 # 20 draws of 10 records: component 2 holds 7 of them, and exchanging the
 # labels lowers the log posterior by about 9 in every draw (a component-2
 # prior that the majority contradicts): component 1 stays the minority
 S <- 20; N <- 10
 z <- matrix(rep(c(1L, 1L, 1L, rep(2L, 7)), each = S), S, N)
 b1 <- matrix(2, S, 2); b2 <- matrix(0, S, 2)
 ex <- rep(-9, S)
 expect_warning(out <- postlink:::align_mixture_labels(z, list(beta = list(b1, b2)), theta = rep(0.3, S),
                                                       orientation = "majority component", relabel = TRUE,
                                                       exchange_lp = ex),
                "estimated from the stored draws, is lower by about 9.0 (more than 2), so the draws were not exchanged after sampling",
                fixed = TRUE)
 expect_true(out$relabelled)
 expect_false(out$flipped)
 expect_true(out$flip_refused)
 expect_equal(out$flip_lp_change, -9)
 expect_false(any(out$swapped))
 expect_identical(out$z, z)
 # within the tolerance (or without exchange changes, i.e. identical priors)
 # the labels are exchanged so that component 1 is the majority
 for (e in list(rep(-1.5, S), NULL)) {
  expect_silent(out <- postlink:::align_mixture_labels(z, list(beta = list(b1, b2)), theta = rep(0.3, S),
                                                       orientation = "majority component", relabel = TRUE,
                                                       exchange_lp = e))
  expect_true(out$flipped)
  expect_false(out$flip_refused)
  expect_true(all(out$swapped))
  expect_equal(out$theta, rep(0.7, S))
  expect_equal(out$pairs$beta[[1L]], b2)
 }
 # the estimate refers to the aligned draws: a draw exchanged by the ECR step
 # contributes the opposite change
 z2 <- z
 z2[1:5, ] <- 3L - z2[1:5, ]
 ex2 <- c(rep(9, 5), rep(-9, S - 5))
 expect_warning(out <- postlink:::align_mixture_labels(z2, list(beta = list(b1, b2)), theta = rep(0.3, S),
                                                       orientation = "majority component", relabel = TRUE,
                                                       exchange_lp = ex2),
                "lower by about 9.0")
 expect_equal(sum(out$swapped), 5L)
 expect_equal(out$flip_lp_change, -9)
 # a refusal before sampling is reported with that of the draws
 m <- postlink:::.minority_message(0.3, list(refused = TRUE, lp_change = -8, flip_refused = TRUE,
                                             flip_lp_change = -9))
 expect_match(m, paste0("lower by about 8.0 (more than 2), so the sampler did not exchange the labels of its ",
                        "starting state, and the draws were not exchanged after sampling either ",
                        "(diagnostics$orient_refused; see 'Label switching' in ?glmMixBayes)."), fixed = TRUE)
 # a refusal after sampling only is not mistaken for one before sampling
 m2 <- postlink:::.minority_message(0.3, list(refused = FALSE, flip_refused = TRUE, flip_lp_change = -9,
                                              default_priors = TRUE))
 expect_match(m2, "estimated from the stored draws, is lower by about 9.0 (more than 2)", fixed = TRUE)
 expect_false(grepl("orient_refused|starting state", m2))
 expect_match(m2, "These are the default priors", fixed = TRUE)
})

test_that("a chain that drifts to the labelling favoured by the component priors is not flipped against them", {
 # lifem records without safe matches and the priors of the example of
 # ?adjMixBayes: beta2 = normal(0, 0.01) says that component 2 has no trend.
 # The flat group holds the majority, so the global exchange after sampling
 # would put the cubic trend into component 2 against its prior (and report
 # a labelling of posterior probability about 0); it is refused
 data(lifem, package = "postlink", envir = environment())
 lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
 d <- rbind(head(subset(lifem, hndlnk == 1), 100), head(subset(lifem, hndlnk == 0), 20))
 pr <- list(intercept1 = "normal(60, 20)", intercept2 = "normal(60, 20)",
            beta1 = "normal(0, 100)", beta2 = "normal(0, 0.01)", theta = "beta(2, 2)")
 q <- with_warnings(plglm(age_at_death ~ poly(unit_yob, 3, raw = TRUE), family = "gaussian",
                          adjustment = adjMixBayes(d, priors = pr),
                          control = list(iterations = 2000, burnin.iterations = 1000, seed = 1)))
 f <- q$value
 dg <- f$diagnostics
 expect_equal(dg$orientation, "majority component")
 expect_true(dg$relabelled)
 # the component-2 draws respect their prior: no trend
 expect_lt(max(abs(colMeans(f$estimates$coefficients2)[-1])), 0.1)
 expect_lt(dg$other_labelling, 0.01)
 if (mean(f$m_samples == 1L) < 0.5) {
  expect_true(any(grepl("favour the labelling in which component 1 is the minority", q$warnings)))
 }
})

test_that("under the match-rate prior the starting state moves to the more probable labelling (regression)", {
 # 65% of the records near 0 and 35% near 4, with m.rate = 0.3 (theta ~
 # beta(14, 6)): the labelling with component 1 on the 35% cluster is about
 # exp(6) times less probable, but the pilot chains of seed 17 all reach it.
 # The main chain used to stay there (theta 0.34, component-1 intercept 4.1)
 # without any warning
 set.seed(4)
 N <- 200
 x <- stats::rnorm(N)
 zt <- ifelse(stats::runif(N) < 0.65, 1, 2)
 y <- ifelse(zt == 1, 0.2 * x + stats::rnorm(N), 4 + stats::rnorm(N))
 X <- cbind("(Intercept)" = 1, x = x)
 ctl <- list(iterations = 2000, burnin.iterations = 1000, seed = 17)
 q <- with_warnings(glmMixBayes(X, y, "gaussian", m.rate = 0.3, control = ctl))
 f <- q$value
 d <- f$diagnostics
 expect_equal(d$orientation, "match-rate prior")
 expect_true(d$pre_oriented)
 expect_gt(d$orient_lp_change, 3)
 expect_false(d$relabelled)
 expect_equal(d$n_swapped, 0L)
 expect_lt(abs(mean(f$estimates$coefficients[, "(Intercept)"])), 0.5)
 expect_gt(mean(f$estimates$theta), 0.6)
 expect_gt(mean((colMeans(f$m_samples == 1L) > 0.5) == (zt == 1)), 0.9)
 expect_lt(d$other_labelling, 0.01)
 expect_length(q$warnings, 0L)
 # the same chain without the orientation stays in the less probable labelling
 ep <- engine_priors(NULL, "gaussian", m.rate = 0.3)
 post <- postlink:::run_mixbayes_engine("gaussian", X, y, NULL, NULL, integer(N), ep,
                                        postlink:::.mixbayes_settings(list(), ctl), pre_orient = "none")
 expect_false(post$pre_oriented)
 expect_lt(mean(post$theta), 0.5)
 expect_gt(postlink:::.other_labelling_share(post$exchange_lp), 0.99)
 # it would be warned about
 expect_warning(postlink:::.check_label_mass(postlink:::.other_labelling_share(post$exchange_lp),
                                             "match-rate prior"),
                "component 1 probably describes the mismatches")

 # a prior close to symmetry (m.rate = 0.45, m.rate.sd = 0.2) identifies the
 # components only weakly, which is reported
 q2 <- with_warnings(glmMixBayes(X, y, "gaussian", m.rate = 0.45, m.rate.sd = 0.2, control = ctl))
 expect_equal(q2$value$diagnostics$orientation, "match-rate prior")
 expect_gt(q2$value$diagnostics$other_labelling, 0.05)
 expect_true(any(grepl("identifies the components only weakly|less probable of the two labellings", q2$warnings)))
})
