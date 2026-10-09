# Parameterisation of the Bayesian mixture models: centring of the linkage
# covariates of the match-probability model (the gamma_intercept prior refers
# to a record with average linkage covariates), the identified Weibull
# intercept (Intercept) + log(scale), and offset() terms in the model formula.
local_edition(3)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Gaussian mixture whose match probability depends on an uncentred linkage
# covariate. `w0` takes the values -1.5, -0.5, 0.5, 1.5 equally often among
# the records that are not flagged as safe matches, so its mean over them is
# exactly 0 in floating point; the safe matches (if any) have w0 = 1.5.
make_link <- function(n = 240, n_safe = 0L, seed = 11) {
 set.seed(seed)
 w0 <- sample(rep(c(-1.5, -0.5, 0.5, 1.5), length.out = n))
 match <- stats::rbinom(n, 1, stats::plogis(1.5 + 1.5 * w0)) == 1
 safe <- rep(0L, n)
 if (n_safe > 0L) {
  w0 <- c(w0, rep(1.5, n_safe))
  match <- c(match, rep(TRUE, n_safe))
  safe <- c(safe, rep(1L, n_safe))
 }
 x <- stats::rnorm(length(w0))
 y <- ifelse(match, 1 + 2 * x, -2 - x) + stats::rnorm(length(w0), sd = 0.5)
 list(X = cbind("(Intercept)" = 1, x = x), y = y, w0 = w0, safe = safe)
}

# Draws of gamma from the engine run on the UNCENTRED linkage covariates,
# aligned as glmMixBayes() aligns them (the pipeline before the covariates
# were centred).
uncentred_gamma <- function(X, y, Z, priors, ctl) {
 safe <- rep(0L, nrow(X))
 pf <- postlink:::prepare_mixbayes_priors(priors, "gaussian", "glm", use_logistic = TRUE)
 ep <- postlink:::build_engine_priors(pf, "gaussian", "glm", K = ncol(X), M = ncol(Z), intercept = TRUE)
 mc <- postlink:::.mixbayes_settings(list(), ctl)
 rule <- postlink:::.label_rule(ep, safe)
 post <- postlink:::run_mixbayes_engine("gaussian", X, y, NULL, Z, safe, ep, mc,
                                        pre_orient = rule$pre_orient)
 al <- postlink:::align_mixture_labels(
  post$z, list(beta = list(post$beta1, post$beta2), disp = list(post$disp1, post$disp2)),
  gamma = post$gamma, safe = safe, orientation = rule$orientation, relabel = rule$exchangeable,
  expect_majority = postlink:::.prior_majority(pf, TRUE))
 al$gamma
}

# Raw (sampled) draws of a Weibull fit: the engine and the label alignment of
# survregMixBayes(), without the identified-intercept transformation.
raw_weibull <- function(X, time, event, ctl, priors = NULL) {
 safe <- rep(0L, nrow(X))
 pf <- postlink:::prepare_mixbayes_priors(priors, "weibull", "survival")
 ep <- postlink:::build_engine_priors(pf, "weibull", "survival", K = ncol(X),
                                      intercept = postlink:::.has_intercept_column(X))
 mc <- postlink:::.mixbayes_settings(list(), ctl)
 rule <- postlink:::.label_rule(ep, safe)
 post <- postlink:::run_mixbayes_engine("surv_weibull", X, time, event, NULL, safe, ep, mc,
                                        pre_orient = rule$pre_orient)
 al <- postlink:::align_mixture_labels(
  post$z, list(beta = list(post$beta1, post$beta2), shape = list(post$disp1, post$disp2),
               scale = list(post$scale1, post$scale2)),
  theta = post$theta, safe = safe, orientation = rule$orientation, relabel = rule$exchangeable,
  expect_majority = postlink:::.prior_majority(pf, FALSE))
 list(b1 = al$pairs$beta[[1L]], b2 = al$pairs$beta[[2L]],
      scale1 = al$pairs$scale[[1L]], scale2 = al$pairs$scale[[2L]])
}

make_weibull <- function(n = 160, seed = 5) {
 set.seed(seed)
 x <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, 0.8) == 1
 tt <- stats::rweibull(n, shape = 1.8, scale = exp(ifelse(match, 0.5 + 0.8 * x, 1.5)))
 list(x = x, time = tt, event = rep(1L, n), match = match)
}

ctl_w <- list(iterations = 400, burnin.iterations = 200, seed = 21)

# ------------------------------------------------------------------------------
# A. Centring of the linkage covariates
# ------------------------------------------------------------------------------

test_that("the centring helpers use the records not flagged as safe and map gamma back", {
 Z <- cbind("(Intercept)" = 1, a = c(1, 2, 3, 10), b = c(0, 4, 8, -100))
 safe <- c(0L, 0L, 0L, 1L)
 zc <- postlink:::.z_center(Z, safe)
 expect_equal(zc, c(a = 2, b = 4))                       # safe record excluded
 expect_equal(postlink:::.z_center(Z, rep(1L, 4)), colMeans(Z[, -1]))  # all safe: all records
 expect_length(postlink:::.z_center(Z[, 1, drop = FALSE], safe), 0L)
 Zc <- postlink:::.center_Z(Z, zc)
 expect_equal(Zc[, 1], rep(1, 4))
 expect_equal(colMeans(Zc[safe == 0, -1]), c(a = 0, b = 0))

 # starting values on the original scale -> centred model, and back
 g <- c(0.3, -1.2, 0.25)
 gc <- postlink:::.center_init_gamma(list(gamma = g, beta1 = 1), zc)
 expect_equal(gc$beta1, 1)
 expect_equal(drop(Zc %*% gc$gamma), drop(Z %*% g))     # same linear predictor
 expect_equal(gc$gamma[-1], g[-1])
 draws <- rbind(gc$gamma, 2 * gc$gamma)
 back <- postlink:::.uncenter_gamma(draws, zc)
 expect_equal(back[, 1], draws[, 1] - drop(draws[, -1] %*% zc))
 expect_equal(back[1, ], g)
 expect_equal(drop(Z %*% t(back)), drop(Zc %*% t(draws)))
 # invalid starting values are left for the validation of control$init
 expect_identical(postlink:::.center_init_gamma(list(gamma = 1), zc), list(gamma = 1))
 expect_identical(postlink:::.center_init_gamma(NULL, zc), NULL)
})

test_that("gamma is reported for the original covariates, with init on the original scale", {
 d <- make_link(n_safe = 20L)
 ZA <- cbind("(Intercept)" = 1, w = d$w0 + 2)   # uncentred: mean 2 over the non-safe records
 ZB <- cbind("(Intercept)" = 1, w = d$w0)       # already centred
 ctl <- list(iterations = 300, burnin.iterations = 150, seed = 8)
 fit <- function(Z, init = NULL) {
  suppressWarnings(glmMixBayes(d$X, d$y, "gaussian", Z = Z, safe.matches = d$safe,
                               control = c(ctl, if (!is.null(init)) list(init = list(gamma = init)))))
 }
 fA <- fit(ZA)
 fB <- fit(ZB)
 # the safe matches (w0 = 1.5) are left out of the centre
 expect_equal(fA$z_center, c(w = 2))
 expect_equal(fB$z_center, c(w = 0))
 expect_equal(unname(fA$z_center), mean(ZA[d$safe == 0, "w"]))
 expect_false(isTRUE(all.equal(unname(fA$z_center), mean(ZA[, "w"]))))
 # the sampler sees the same centred covariates: identical draws, mapped back
 # to each Z by intercept = intercept_c - slope * centre
 expect_equal(fA$estimates$coefficients, fB$estimates$coefficients)
 expect_equal(fA$estimates$gamma[, "w"], fB$estimates$gamma[, "w"])
 expect_equal(fA$estimates$gamma[, "(Intercept)"],
              fB$estimates$gamma[, "(Intercept)"] - 2 * fB$estimates$gamma[, "w"])
 expect_equal(fA$estimates$m.coefficients, -fA$estimates$gamma)
 expect_equal(colnames(fA$estimates$gamma), c("(Intercept)", "w"))

 # starting values are given for the original covariates: c(a, b) for ZA is
 # c(a + 2 b, b) for ZB, so the two chains coincide
 iA <- fit(ZA, init = c(-1, 0.8))
 iB <- fit(ZB, init = c(-1 + 2 * 0.8, 0.8))
 expect_equal(iA$diagnostics$settings$pilots, 0L)
 expect_equal(iA$estimates$coefficients, iB$estimates$coefficients)
 expect_equal(iA$estimates$gamma[, "(Intercept)"],
              iB$estimates$gamma[, "(Intercept)"] - 2 * iB$estimates$gamma[, "w"])
 expect_error(fit(ZA, init = c(1, 2, 3)), "control\\$init\\$gamma")
})

test_that("with a weak intercept prior the posterior of gamma does not depend on the centring", {
 d <- make_link(n = 400, seed = 12)
 Z <- cbind("(Intercept)" = 1, w = d$w0 + 3)
 pri <- list(gamma_intercept = "normal(0, 100)")
 ctl <- list(iterations = 3000, burnin.iterations = 1000, seed = 3)
 fc <- suppressWarnings(glmMixBayes(d$X, d$y, "gaussian", Z = Z, priors = pri, control = ctl))
 gu <- suppressWarnings(uncentred_gamma(d$X, d$y, Z, pri, ctl))
 gc <- fc$estimates$gamma
 mcse <- function(v) stats::sd(v) / sqrt(postlink:::.ess_ar(v))
 for (j in 1:2) {
  tol <- 4 * sqrt(mcse(gc[, j])^2 + mcse(gu[, j])^2)
  expect_lt(abs(mean(gc[, j]) - mean(gu[, j])), tol, label = colnames(gc)[j])
 }
 # and the slope is recovered (true value 1.5)
 expect_lt(abs(mean(gc[, "w"]) - 1.5), 0.6)
})

test_that("the prior implied by m.rate refers to the average record, as print.adjMixBayes() says", {
 set.seed(2)
 n <- 150
 dat <- data.frame(w = stats::rnorm(n, 3), x = stats::rnorm(n))
 dat$y <- 1 + dat$x + stats::rnorm(n, sd = 0.5)
 dat$safe <- TRUE        # every record is safe: the match model keeps its prior
 adj <- adjMixBayes(dat, m.formula = ~ w, m.rate = 0.1, m.rate.sd = 0.03, safe.matches = safe)
 out <- capture_output(print(adj))
 bp <- postlink:::mrate_to_beta(0.1, 0.03)
 mom <- postlink:::logit_beta_moments(bp$alpha, bp$beta)
 q <- stats::plogis(mom[["mu"]] + c(0, -1, 1) * stats::qnorm(0.975) * mom[["sd"]])
 expect_match(out, sprintf("Implied Prior on theta: median %s, 95%% interval \\(%s, %s\\)",
                           format(q[1], digits = 3), format(q[2], digits = 3), format(q[3], digits = 3)))
 expect_match(out, "for a record with average linkage covariates")
 expect_match(out, sprintf("w = %s", format(mean(dat$w), digits = 3)), fixed = TRUE)
 expect_match(out, "the intercept refers to a record with average linkage covariates")
 # every record is safe: the mean is over all records, and the print says so
 expect_match(out, "average linkage covariates over all records
      (every record is flagged as a safe match)",
              fixed = TRUE)
 expect_match(out, "centred at their mean over all records, as every record is flagged", fixed = TRUE)
 expect_false(grepl("not flagged as safe", out, fixed = TRUE))
 expect_true(attr(postlink:::.adj_z_center(adj), "all_safe"))

 fit <- suppressWarnings(plglm(y ~ x, family = "gaussian", adjustment = adj,
                               control = list(iterations = 4000, burnin.iterations = 1000, seed = 5)))
 expect_equal(unname(fit$z_center), mean(dat$w))
 g <- fit$estimates$gamma
 # theta of a record with average linkage covariates
 th <- stats::plogis(g[, "(Intercept)"] + g[, "w"] * fit$z_center)
 expect_equal(unname(stats::quantile(th, c(0.5, 0.025, 0.975))), q, tolerance = 0.02)
 # the gamma_intercept prior derived from m.rate is that of the centred model
 expect_equal(unname(fit$priors["gamma_intercept"]),
              sprintf("normal(%s, %s)", format(signif(mom[["mu"]], 4)), format(signif(mom[["sd"]], 4))))

 # without linked data the mean is computed at fit time
 out2 <- capture_output(print(adjMixBayes(NULL, m.formula = ~ w, m.rate = 0.1, m.rate.sd = 0.03)))
 expect_match(out2, "computed when the model is fitted")
})

test_that("print.adjMixBayes() uses the records not flagged as safe for the average record", {
 dat <- data.frame(w = c(1, 2, 3, 100), y = 1:4)
 adj <- adjMixBayes(dat, m.formula = ~ w, m.rate = 0.2, safe.matches = c(FALSE, FALSE, FALSE, TRUE))
 expect_equal(postlink:::.adj_z_center(adj), structure(c(w = 2), all_safe = FALSE))
 out <- capture_output(print(adj))
 expect_match(out, "average linkage covariates over the records not flagged as safe matches", fixed = TRUE)
 expect_match(out, "(in linked.data: w = 2; fit$z_center holds the mean over the analysed records)",
              fixed = TRUE)
 expect_false(grepl("every record is flagged", out, fixed = TRUE))
 expect_null(postlink:::.adj_z_center(adjMixBayes(dat)))
})

test_that("summary() states that the gamma_intercept prior applies at z_center", {
 set.seed(7)
 n <- 100
 d <- data.frame(x = stats::rnorm(n), w = stats::rnorm(n, 5), w2 = stats::runif(n))
 d$y <- 1 + d$x + stats::rnorm(n, sd = 0.5)
 d$time <- stats::rweibull(n, 2, exp(0.3 + 0.4 * d$x))
 d$status <- 1L
 ctl <- list(iterations = 100, burnin.iterations = 50, seed = 3)
 adj <- adjMixBayes(d, m.formula = ~ w + w2, m.rate = 0.2)
 note <- function(f) {
  sprintf("(the gamma_intercept prior %s applies at the average linkage covariates
 z_center: w = %s, w2 = %s;",
          f$priors[["gamma_intercept"]], format(f$z_center[["w"]], digits = 4),
          format(f$z_center[["w2"]], digits = 4))
 }
 fg <- suppressWarnings(plglm(y ~ x, adjustment = adj, control = ctl))
 sg <- summary(fg)
 expect_equal(sg$z_center, fg$z_center)
 expect_identical(sg$gamma_intercept_prior, unname(fg$priors[["gamma_intercept"]]))
 expect_match(paste(utils::capture.output(print(sg, digits = 4)), collapse = "
"), note(fg), fixed = TRUE)
 fs <- suppressWarnings(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull", adjustment = adj,
                                  control = ctl))
 ss <- summary(fs)
 expect_equal(ss$z_center, fs$z_center)
 expect_match(paste(utils::capture.output(print(ss, digits = 4)), collapse = "
"), note(fs), fixed = TRUE)

 # no note without linkage covariates, nor for fits saved without z_center
 f0 <- suppressWarnings(plglm(y ~ x, adjustment = adjMixBayes(d), control = ctl))
 expect_null(summary(f0)$z_center)
 expect_false(grepl("z_center", paste(utils::capture.output(print(summary(f0))), collapse = "
")))
 fg$z_center <- NULL
 expect_false(grepl("z_center", paste(utils::capture.output(print(summary(fg))), collapse = "
")))
 # unnamed linkage covariates are labelled by their column of Z
 expect_output(postlink:::.print_z_center_note(c(0.5, 2), NULL, 3),
               "applies at the average linkage covariates
 z_center: Z[, 2] = 0.5, Z[, 3] = 2;", fixed = TRUE)
})

test_that("survregMixBayes() centres the linkage covariates too", {
 d <- make_weibull()
 X <- cbind("(Intercept)" = 1, x = d$x)
 set.seed(9)
 w <- stats::rnorm(length(d$x), 5)
 Z <- cbind("(Intercept)" = 1, w = w)
 safe <- rep(0L, length(w)); safe[1:10] <- 1L
 f <- suppressWarnings(survregMixBayes(X, cbind(d$time, d$event), dist = "weibull", Z = Z,
                                       safe.matches = safe, control = ctl_w))
 expect_equal(f$z_center, c(w = mean(w[safe == 0])))
 expect_equal(f$estimates$m.coefficients, -f$estimates$gamma)
 # the reported gamma gives the same linear predictor as the centred model
 gcen <- cbind(f$estimates$gamma[, 1] + f$estimates$gamma[, 2] * f$z_center, f$estimates$gamma[, 2])
 expect_equal(Z %*% t(f$estimates$gamma), postlink:::.center_Z(Z, f$z_center) %*% t(gcen))
})

# ------------------------------------------------------------------------------
# B. Identified Weibull intercept
# ------------------------------------------------------------------------------

test_that("Weibull fits report (Intercept) + log(scale) draw by draw", {
 d <- make_weibull()
 X <- cbind("(Intercept)" = 1, x = d$x)
 f <- suppressWarnings(survregMixBayes(X, cbind(d$time, d$event), dist = "weibull", control = ctl_w))
 r <- suppressWarnings(raw_weibull(X, d$time, d$event, ctl_w))
 expect_true(f$intercept_includes_logscale)
 b <- f$estimates$coefficients
 expect_equal(unname(b[, "(Intercept)"]), unname(r$b1[, 1] + log(r$scale1)))
 expect_equal(unname(b[, "x"]), unname(r$b1[, 2]))
 expect_equal(unname(f$estimates$coefficients2[, 1]), unname(r$b2[, 1] + log(r$scale2)))
 expect_equal(unname(f$estimates$coefficients2[, 2]), unname(r$b2[, 2]))
 expect_equal(f$estimates$scale, r$scale1)
 expect_equal(f$estimates$scale2, r$scale2)

 # predictions are those of the raw parameterisation (log(scale) added once)
 nx <- X[1:5, ]
 expect_equal(unname(predict(f, newdata = nx)$component1),
              rowMeans(nx %*% t(r$b1) + matrix(log(r$scale1), 5, nrow(r$b1), byrow = TRUE)))
 expect_equal(unname(predict(f, newdata = nx)$component2),
              rowMeans(nx %*% t(r$b2) + matrix(log(r$scale2), 5, nrow(r$b2), byrow = TRUE)))

 # every method uses the identified intercept
 expect_equal(coef(f), colMeans(b))
 expect_equal(vcov(f), stats::cov(b))
 expect_equal(confint(f, block = "coefficients")[, 1],
              apply(b, 2, stats::quantile, probs = 0.025, names = FALSE))
 expect_equal(unname(posterior_draws(f)[, "coefficients[(Intercept)]"]), unname(b[, "(Intercept)"]))
 expect_equal(unname(f$diagnostics$ess["coefficients[(Intercept)]"]),
              postlink:::.ess_ar(b[, "(Intercept)"]))
 s <- summary(f)
 expect_equal(unname(s$coef1["(Intercept)", "Estimate"]), mean(b[, "(Intercept)"]))
 expect_null(s$intercept_logscale)
 expect_equal(s$inv_shape1[, 1:4, drop = FALSE], postlink:::.posterior_table(1 / f$estimates$shape, "1/shape"))
 expect_equal(unname(s$inv_shape1[, "Rhat"]), postlink:::.rhat_split(1 / f$estimates$shape))
 out <- paste(utils::capture.output(print(s)), collapse = "\n")
 expect_match(out, "comparable to the Scale of survival::survreg()", fixed = TRUE)
 expect_match(out, "not identified separately from the intercept: reflects its prior", fixed = TRUE)
 expect_match(out, "the identified (Intercept) + log(scale)", fixed = TRUE)
 expect_match(paste(utils::capture.output(print(f)), collapse = "\n"),
              "identified Weibull intercept", fixed = TRUE)
})

test_that("the intercept column is found by its values, and designs without one are unchanged", {
 d <- make_weibull()
 # intercept in the second column, without names
 X2 <- cbind(d$x, 1)
 f2 <- suppressWarnings(survregMixBayes(X2, cbind(d$time, d$event), dist = "weibull", control = ctl_w))
 r2 <- suppressWarnings(raw_weibull(X2, d$time, d$event, ctl_w))
 expect_true(f2$intercept_includes_logscale)
 expect_equal(unname(f2$estimates$coefficients[, 2]), unname(r2$b1[, 2] + log(r2$scale1)))
 expect_equal(unname(f2$estimates$coefficients[, 1]), unname(r2$b1[, 1]))

 # no intercept column: reported as sampled, predict() adds log(scale)
 X0 <- cbind(x = d$x)
 f0 <- suppressWarnings(survregMixBayes(X0, cbind(d$time, d$event), dist = "weibull", control = ctl_w))
 r0 <- suppressWarnings(raw_weibull(X0, d$time, d$event, ctl_w))
 expect_false(f0$intercept_includes_logscale)
 expect_equal(unname(f0$estimates$coefficients), unname(r0$b1))
 expect_equal(unname(f0$estimates$coefficients2), unname(r0$b2))
 nx <- X0[1:4, , drop = FALSE]
 expect_equal(unname(predict(f0, newdata = nx)$component1),
              rowMeans(nx %*% t(r0$b1) + matrix(log(r0$scale1), 4, nrow(r0$b1), byrow = TRUE)))
 s0 <- summary(f0)
 expect_false(s0$intercept_includes_logscale)
 out0 <- paste(utils::capture.output(print(s0)), collapse = "\n")
 expect_false(grepl("not identified separately", out0, fixed = TRUE))
 expect_match(out0, "Scale multiplier of exp(X beta) (component 1 = correct-match)", fixed = TRUE)

 # gamma survival fits have no scale and no such element
 fg <- suppressWarnings(survregMixBayes(cbind(1, d$x), cbind(d$time, d$event), dist = "gamma", control = ctl_w))
 expect_null(fg$intercept_includes_logscale)
 expect_null(summary(fg)$inv_shape1)
})

test_that("the prior-data conflict check uses the sampled Weibull intercept", {
 d <- make_weibull(seed = 6)
 X <- cbind("(Intercept)" = 1, x = d$x)
 tt <- d$time * exp(2)      # identified intercept about 2.5
 # safe matches anchor component 1 (no relabelling)
 safe <- as.integer(d$match & seq_along(d$match) <= 100)
 catch <- function(expr) {
  msgs <- character(0)
  res <- withCallingHandlers(expr, warning = function(w) {
   msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  list(fit = res, warnings = msgs)
 }
 # a tight intercept prior at 0: the sampled intercept stays near 0 (log(scale)
 # absorbs the level), so there is no conflict, although the identified
 # intercept lies far from 0
 r <- catch(survregMixBayes(X, cbind(tt, d$event), dist = "weibull", safe.matches = safe,
                            priors = list(intercept1 = "normal(0, 0.05)"), control = ctl_w))
 expect_gt(mean(r$fit$estimates$coefficients[, "(Intercept)"]), 1.5)
 expect_false(any(grepl("prior standard deviations from its prior mean", r$warnings)))
 # with the scale fixed near 1 the sampled intercept must move: conflict
 r2 <- catch(survregMixBayes(X, cbind(tt, d$event), dist = "weibull", safe.matches = safe,
                             priors = list(intercept1 = "normal(0, 0.2)", scale1 = "lognormal(0, 0.01)"),
                             control = ctl_w))
 expect_true(any(grepl("'\\(Intercept\\)' in the correct-match component", r2$warnings)))
})

# ------------------------------------------------------------------------------
# C. offset() terms in the model formula
# ------------------------------------------------------------------------------

test_that("offset() terms are refused for adjMixBayes and left alone for the frequentist fits", {
 set.seed(4)
 n <- 80
 d <- data.frame(x = stats::rnorm(n), e = stats::runif(n, 1, 3), x2 = stats::rnorm(n))
 d$y <- stats::rpois(n, d$e * exp(0.2 + 0.5 * d$x))
 d$time <- stats::rweibull(n, 2, exp(0.5 + 0.3 * d$x))
 d$status <- 1L
 msg <- "offset() terms in the model formula are not supported by the Bayesian mixture models (adjMixBayes)."
 expect_error(plglm(y ~ x + offset(log(e)), family = "poisson", adjustment = adjMixBayes(d),
                    control = list(iterations = 50)), msg, fixed = TRUE)
 expect_error(plsurvreg(survival::Surv(time, status) ~ x + offset(x2), adjustment = adjMixBayes(d),
                        control = list(iterations = 50)), msg, fixed = TRUE)
 # the frequentist mixture fit is not affected by the check
 expect_s3_class(suppressWarnings(plglm(y ~ x + offset(log(e)), family = "poisson",
                                        adjustment = adjMixture(d))), "glmMixture")
})

test_that("unsupported modelling arguments are refused by name, before they are evaluated", {
 set.seed(4)
 n <- 60
 d <- data.frame(x = stats::rnorm(n), e = stats::runif(n, 1, 3))
 d$y <- stats::rpois(n, d$e * exp(0.2 + 0.5 * d$x))
 d$time <- stats::rweibull(n, 2, exp(0.5 + 0.3 * d$x))
 d$status <- 1L
 ctl <- list(iterations = 50)
 # `e` and `not_defined` exist only as columns of the linked data / nowhere:
 # the arguments cannot be evaluated, but they are refused before that
 expect_false(exists("not_defined"))
 msg <- "The Bayesian mixture models do not support the argument(s) offset."
 expect_error(plglm(y ~ x, family = "poisson", offset = log(e), adjustment = adjMixBayes(d), control = ctl),
              msg, fixed = TRUE)
 expect_error(plsurvreg(survival::Surv(time, status) ~ x, offset = not_defined, adjustment = adjMixBayes(d),
                        control = ctl), msg, fixed = TRUE)
 expect_error(plglm(y ~ x, family = "poisson", weights = e, adjustment = adjMixBayes(d), control = ctl),
              "do not support the argument(s) weights.", fixed = TRUE)
 X <- cbind(1, d$x)
 expect_error(glmMixBayes(X, d$y, "poisson", offset = not_defined, control = ctl), msg, fixed = TRUE)
 expect_error(survregMixBayes(X, cbind(d$time, d$status), "weibull", start = not_defined, control = ctl),
              "do not support the argument(s) start.", fixed = TRUE)
 # names of `...` are read without evaluating the arguments
 expect_identical(postlink:::.dots_names(a = stop("evaluated"), 1, b = not_defined), c("a", "", "b"))
 expect_identical(postlink:::.dots_names(), character(0))
})
