# Regression tests for behaviour fixed after the review of the C++ Gibbs port:
# row alignment in mi_with(), link functions, response and survival-response
# validation, starting values, MCMC settings, prior bookkeeping, label
# orientation and the survMixBayes / glmMixBayes methods.
local_edition(3)

make_rev <- function(n = 160, seed = 11) {
 set.seed(seed)
 x <- stats::rnorm(n)
 z <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, stats::plogis(1.2 + 1.5 * z)) == 1
 y <- ifelse(match, 1 + 2 * x + stats::rnorm(n, sd = 0.7), stats::rnorm(n, sd = 2))
 tt <- stats::rweibull(n, shape = 2, scale = exp(ifelse(match, 0.5 + 0.8 * x, 0.2)))
 cens <- stats::rexp(n, 0.1)
 data.frame(y = y, x = x, z = z, safe = match & seq_len(n) <= 40,
            time = pmin(tt, cens), status = as.integer(tt <= cens))
}

ctl <- list(iterations = 400, burnin.iterations = 200, seed = 5)

test_that("mi_with() aligns `data` with the records used in the fit", {
 d <- make_rev()
 d$y[c(3, 17, 90)] <- NA
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 expect_identical(colnames(fit$m_samples), fit$obs_names)

 p_full <- mi_with(fit, data = d, formula = y ~ x)
 p_cc <- mi_with(fit, data = d[!is.na(d$y), ], formula = y ~ x)
 expect_equal(p_full$coef, p_cc$coef)
 p_reset <- mi_with(fit, data = `rownames<-`(d[!is.na(d$y), ], NULL), formula = y ~ x)
 expect_equal(p_reset$coef, p_cc$coef)
 expect_error(mi_with(fit, data = d[1:100, ], formula = y ~ x), "rows")

 # direct engine fits: `data` must hold one row per analysed record
 keep <- !is.na(d$y)
 f2 <- suppressMessages(glmMixBayes(cbind(1, d$x[keep]), d$y[keep], "gaussian", control = ctl))
 expect_error(mi_with(f2, data = d, formula = y ~ x), "rows")

 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                  adjustment = adjMixBayes(d, m.formula = ~ z, safe.matches = safe),
                                  control = ctl))
 expect_s3_class(mi_with(fS, data = d), "mi_link_pool_survreg")
})

test_that("mi_with() refits gamma models with the log link and reports skipped draws", {
 d <- make_rev()
 d$ypos <- exp(d$y / 3)
 fit <- suppressMessages(plglm(ypos ~ x, family = Gamma(link = "log"), adjustment = adjMixBayes(d),
                               control = ctl))
 pool <- mi_with(fit, data = d, formula = ypos ~ x)
 expect_equal(pool$m, nrow(fit$m_samples))
 # pooled coefficients are on the log scale of the mixture model
 expect_lt(max(abs(pool$coef - colMeans(fit$estimates$coefficients))), 0.3)
 counts <- rowSums(fit$m_samples == 1L)
 expect_warning(mi_with(fit, data = d, formula = ypos ~ x,
                        min_n = as.integer(stats::median(counts)) + 1L),
                "posterior draws were not used")
 # when no draw can be used, the error says why
 expect_error(mi_with(fit, data = d, formula = ypos ~ x, min_n = nrow(d) + 1L),
              "fewer than `min_n`")
 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                  adjustment = adjMixBayes(d), control = ctl))
 expect_error(mi_with(fS, data = d, formula = survival::Surv(time, status) ~ x, dist = "weibull"),
              "First error")
})

test_that("family links are checked and binomial responses may be logical or factors", {
 d <- make_rev()
 d$ypos <- exp(d$y / 3)
 expect_error(plglm(I(y > 0) ~ x, family = binomial(link = "probit"), adjustment = adjMixBayes(d),
                    control = ctl), "logit link")
 expect_error(plglm(ypos ~ x, family = gaussian(link = "log"), adjustment = adjMixBayes(d),
                    control = ctl), "identity link")
 expect_warning(suppressMessages(plglm(ypos ~ x, family = Gamma(), adjustment = adjMixBayes(d),
                                       control = ctl)), "log link")
 fg <- suppressMessages(plglm(ypos ~ x, family = Gamma(link = "log"), adjustment = adjMixBayes(d),
                              control = ctl))
 expect_equal(fg$family, "gamma")

 fl <- suppressMessages(plglm(I(y > 0) ~ x, family = "binomial", adjustment = adjMixBayes(d),
                              control = ctl))
 expect_equal(fl$family, "binomial")
 d$yf <- factor(ifelse(d$y > 0, "yes", "no"))
 ff <- suppressMessages(plglm(yf ~ x, family = binomial, adjustment = adjMixBayes(d), control = ctl))
 expect_equal(unname(ff$estimates$coefficients), unname(fl$estimates$coefficients))

 X <- cbind(1, d$x)
 expect_error(glmMixBayes(X, letters[seq_len(nrow(X))], "gaussian", control = ctl), "numeric")
 expect_error(glmMixBayes(X, cbind(d$y > 0, d$y <= 0), "binomial", control = ctl), "two-column")
 expect_error(glmMixBayes(X, d$yf, "gaussian", control = ctl), "only supported for family = 'binomial'")
 expect_error(glmMixBayes(X, d$y, family = c("gaussian", "poisson"), control = ctl), "single character")
})

test_that("non-finite inputs and malformed Z are rejected", {
 d <- make_rev()
 X <- cbind(1, d$x)
 X2 <- X
 X2[1, 2] <- Inf
 expect_error(glmMixBayes(X2, d$y, "gaussian", control = ctl), "missing or infinite")
 expect_error(survregMixBayes(X2, survival::Surv(d$time, d$status), "weibull", control = ctl),
              "missing or infinite")
 y2 <- d$y
 y2[3] <- Inf
 expect_error(glmMixBayes(X, y2, "gaussian", control = ctl), "missing or infinite")
 Z <- cbind(1, d$z)
 Z2 <- Z
 Z2[2, 2] <- NA
 expect_error(glmMixBayes(X, d$y, "gaussian", Z = Z2, control = ctl), "missing or infinite")
 expect_error(glmMixBayes(X, d$y, "gaussian", Z = Z[, 2, drop = FALSE], control = ctl), "column of ones")
 d$x0 <- abs(d$x)
 d$x0[4] <- 0
 expect_error(plglm(y ~ log(x0), family = "gaussian", adjustment = adjMixBayes(d), control = ctl),
              "infinite")
})

test_that("only right-censored survival responses with 0/1 events are accepted", {
 d <- make_rev()
 X <- cbind(1, d$x)
 expect_error(survregMixBayes(X, survival::Surv(d$time, d$status, type = "left"), "weibull", control = ctl),
              "right-censored")
 expect_error(survregMixBayes(X, survival::Surv(d$time, d$time + 1, type = "interval2"), "weibull",
                              control = ctl), "right-censored")
 expect_error(survregMixBayes(X, survival::Surv(rep(0.01, nrow(d)), d$time + 1, d$status), "weibull",
                              control = ctl), "right-censored")
 expect_error(survregMixBayes(X, cbind(d$time, d$status, 1), "weibull", control = ctl), "2-column")
 expect_error(survregMixBayes(X, cbind(d$time, d$status + 1), "weibull", control = ctl), "coded 0")
 f <- suppressMessages(survregMixBayes(X, cbind(d$time, d$status == 1), " Weibull ", control = ctl))
 expect_equal(f$dist, "weibull")
 d$start <- 0.001
 expect_error(plsurvreg(survival::Surv(start, time, status) ~ x, dist = "weibull",
                        adjustment = adjMixBayes(d), control = ctl), "right-censored")
})

test_that("starting values are validated per model and the Weibull scale can be initialised", {
 d <- make_rev()
 X <- cbind("(Intercept)" = 1, x = d$x)
 S <- survival::Surv(d$time, d$status)
 f <- suppressMessages(survregMixBayes(X, S, "weibull", control = c(ctl, list(init = list(
  beta1 = c(0, 0.5), beta2 = c(0, 0), disp1 = 2, disp2 = 1.5, scale1 = 1, scale2 = 1, theta = 0.8)))))
 expect_s3_class(f, "survMixBayes")
 expect_equal(f$diagnostics$settings$pilots, 0L)
 expect_error(survregMixBayes(X, S, "weibull", control = c(ctl, list(init = list(scale1 = -1)))), "positive")
 expect_error(survregMixBayes(X, S, "gamma", control = c(ctl, list(init = list(scale1 = 1)))), "weibull")
 expect_error(glmMixBayes(X, d$y, "gaussian", control = c(ctl, list(init = list(theta = 1.5)))),
              "between 0 and 1")
 expect_error(glmMixBayes(X, d$y, "gaussian", control = c(ctl, list(init = list(disp1 = 0)))), "positive")
 expect_error(glmMixBayes(X, d$y, "gaussian", control = c(ctl, list(init = list(gamma = 1)))), "use `theta`")
 expect_error(glmMixBayes(X, as.integer(d$y > 0), "binomial", control = c(ctl, list(init = list(disp1 = 1)))),
              "dispersion")
 expect_error(glmMixBayes(X, d$y, "gaussian", Z = cbind(1, d$z),
                          control = c(ctl, list(init = list(theta = 0.5)))), "use `gamma`")
 fi <- suppressMessages(glmMixBayes(X, d$y, "gaussian", control = c(ctl, list(init = list(theta = 0.7)))))
 expect_true(all(is.finite(fi$estimates$theta)))
})

test_that("MCMC settings and unsupported modelling arguments are validated", {
 d <- make_rev()
 X <- cbind(1, d$x)
 expect_error(glmMixBayes(X, d$y, "gaussian", control = list(100, 50)), "named list")
 expect_error(glmMixBayes(X, d$y, "gaussian", control = c(ctl, list(collapse = "yes"))), "TRUE or FALSE")
 expect_error(glmMixBayes(X, d$y, "gaussian", control = list(iterations = 400, burnin.iterations = 200,
                                                             seed = "a")), "single integer")
 expect_error(glmMixBayes(X, d$y, "gaussian", control = list(iterations = Inf, burnin.iterations = 200)),
              "single integer")
 expect_warning(suppressMessages(glmMixBayes(X, d$y, "gaussian",
                                             control = c(ctl, list(adapt_delta = 0.9, cores = 2)))),
                "adapt_delta")
 f0 <- suppressMessages(glmMixBayes(X, d$y, "gaussian", control = c(ctl, list(pilots = 0, pilot.iterations = 0))))
 expect_equal(f0$diagnostics$settings$pilots, 0L)
 expect_type(f0$diagnostics$settings$iterations, "integer")
 expect_error(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl,
                    weights = rep(1, nrow(d))), "do not support")
 expect_warning(fd <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d),
                                             control = ctl, data = d)),
                "is ignored")
 expect_s3_class(fd, "glmMixBayes")
 expect_error(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), subset = x > 100,
                    control = ctl), "No observations")
 # 0.1.2 argument order: control is the fifth argument
 f5 <- suppressMessages(glmMixBayes(X, d$y, "gaussian", NULL, ctl))
 expect_equal(f5$diagnostics$settings$seed, 5L)
})

test_that("priors are validated early, stored in the fit and reported", {
 d <- make_rev()
 expect_error(adjMixBayes(d, priors = list(theta = "uniform(0,1)")), "not supported")
 expect_error(adjMixBayes(d, priors = list(beta1 = "normal(0)")), "expected 2 arguments")
 expect_warning(adjMixBayes(d, priors = list(beta3 = "normal(0,1)")), "unrecognised")
 expect_warning(adjMixBayes(d, m.rate.sd = 0.05), "no effect")
 expect_error(adjMixBayes(d, m.rate = 0.01), "strictly less than")

 fit <- suppressMessages(plglm(y ~ x, family = "gaussian",
                               adjustment = adjMixBayes(d, m.rate = 0.2, m.rate.sd = 0.05), control = ctl))
 bp <- postlink:::mrate_to_beta(0.2, 0.05)
 expect_equal(fit$priors[["theta"]],
              sprintf("beta(%s, %s)", format(signif(bp$alpha, 4)), format(signif(bp$beta, 4))))
 expect_equal(fit$priors[["sigma1"]], "cauchy(0, 2.5)")
 expect_message(suppressWarnings(plglm(y ~ x, family = "gaussian",
                                       adjustment = adjMixBayes(d, m.rate = 0.2, priors = list(theta = "beta(2,2)")),
                                       control = ctl)), "not used")
 # (the beta prior of m.rate = 0.05 also identifies the components only
 # weakly on these data, which a second warning reports)
 ws <- character()
 withCallingHandlers(suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d, m.rate = 0.05),
                                            control = ctl)),
                     warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
 expect_true(any(grepl("both shape parameters exceed 1", ws)))

 # priors given at fit time replace those of the adjustment object
 adjS <- adjMixBayes(d, safe.matches = safe, priors = list(intercept1 = "normal(-50, 0.01)"))
 # (the fixed intercept 50 keeps every other record out of component 1, hence
 # the warning; the tight prior also contradicts the safe matches)
 expect_warning(
  expect_warning(fB <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjS, control = ctl,
                                              priors = list(intercept1 = "normal(50, 0.01)"))),
                 "safe matches may be atypical"),
  "The supplied prior of '(Intercept)' in the correct-match component, normal(50, 0.01), conflicts with the data",
  fixed = TRUE)
 expect_equal(fB$priors[["intercept1"]], "normal(50, 0.01)")
 expect_equal(mean(fB$estimates$coefficients[, "(Intercept)"]), 50, tolerance = 0.01)

 # a leading column of ones is the intercept; otherwise all columns get the slope prior
 flat <- postlink:::prepare_mixbayes_priors(NULL, "gaussian", "glm")
 expect_equal(postlink:::build_engine_priors(flat, "gaussian", "glm", K = 2L, intercept = FALSE)$beta1_sd, c(5, 5))
 expect_true(postlink:::.has_intercept_column(cbind(1, 2:3)))
 expect_false(postlink:::.has_intercept_column(cbind(2:3, 1)))
 fN <- suppressMessages(plglm(y ~ 0 + x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 expect_equal(colnames(fN$estimates$coefficients), "x")

 expect_equal(adjMixBayes(d, safe.matches = "safe")$safe.matches, d$safe)
})

test_that("label orientation follows safe matches, the match-rate prior or the majority", {
 set.seed(9)
 n <- 200
 x <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, 0.3) == 1
 y <- ifelse(match, 1 + 2 * x + stats::rnorm(n, sd = 0.5), -3 - x + stats::rnorm(n, sd = 0.5))
 X <- cbind("(Intercept)" = 1, x = x)
 agree <- function(f) mean((colMeans(f$m_samples == 1L) > 0.5) == match)

 f_flat <- suppressMessages(glmMixBayes(X, y, "gaussian", control = ctl))
 expect_equal(f_flat$diagnostics$orientation, "majority component")
 expect_true(f_flat$diagnostics$relabelled)
 expect_gt(mean(f_flat$estimates$theta), 0.5)

 f_prior <- suppressMessages(glmMixBayes(X, y, "gaussian", priors = list(theta = "beta(3,7)"), control = ctl))
 expect_equal(f_prior$diagnostics$orientation, "match-rate prior")
 expect_false(f_prior$diagnostics$relabelled)
 expect_identical(f_prior$diagnostics$pre_oriented, isTRUE(f_prior$diagnostics$orient_lp_change > 0))
 expect_equal(f_prior$diagnostics$n_swapped, 0L)
 expect_lt(f_prior$diagnostics$other_labelling, 0.05)
 expect_lt(mean(f_prior$estimates$theta), 0.45)
 expect_gt(agree(f_prior), 0.9)

 # 70% of the records are mismatches but the default prior does not say so:
 # the safe matches identify the components and a warning points to m.rate
 expect_warning(f_safe <- suppressMessages(glmMixBayes(X, y, "gaussian",
                                                       safe.matches = as.integer(match & seq_len(n) %% 5 == 0),
                                                       control = ctl)),
                "m.rate")
 expect_no_warning(suppressMessages(glmMixBayes(X, y, "gaussian",
                                                safe.matches = as.integer(match & seq_len(n) %% 5 == 0),
                                                m.rate = 0.7, m.rate.sd = 0.1, control = ctl)))
 expect_equal(f_safe$diagnostics$orientation, "safe matches")
 expect_equal(f_safe$diagnostics$n_safe, sum(match & seq_len(n) %% 5 == 0))
 expect_gt(agree(f_safe), 0.9)
 expect_equal(f_safe$diagnostics$n_swapped, 0L)

 # Under the majority convention the sampler orients its starting state so
 # that component 1 holds the majority, unless the component-specific priors
 # make that labelling much less probable: here beta2 says that component 2
 # has no slope, so the flat 70% group is component 2, the slope-2 group stays
 # in component 1, the draws are not exchanged after sampling either (the
 # same priors refuse the global exchange of the aligned draws), and a
 # warning explains why component 1 is the minority
 y0 <- ifelse(match, 1 + 2 * x + stats::rnorm(n, sd = 0.5), 0.5 + stats::rnorm(n, sd = 0.5))
 for (sd2 in c("0.01", "0.1", "0.5")) {
  for (seed in c(1, 5)) {
   ws <- character()
   f_asym <- withCallingHandlers(
    suppressMessages(glmMixBayes(X, y0, "gaussian", priors = list(beta2 = sprintf("normal(0, %s)", sd2)),
                                 control = list(iterations = 1000, burnin.iterations = 500, seed = seed))),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
   d <- f_asym$diagnostics
   expect_equal(d$orientation, "majority component")
   expect_true(d$relabelled)
   expect_equal(d$n_swapped, 0L)
   expect_equal(f_asym$priors[["beta2"]], sprintf("normal(0, %s)", sd2))
   # the main mode: the slope-2 group is recovered in component 1, at the log
   # posterior level of the selected pilot chain
   expect_gt(mean(d$lp), max(d$pilot_lp) - 5)
   expect_lt(abs(mean(f_asym$estimates$coefficients[, "x"]) - 2), 0.2)
   expect_gt(agree(f_asym), 0.8)
   expect_lt(abs(mean(f_asym$estimates$coefficients2[, "x"])), 0.1)
   # the starting state had component 2 as the majority and was not oriented
   expect_true(d$orient_refused)
   expect_false(d$pre_oriented)
   expect_lt(d$orient_lp_change, -postlink:::.ORIENT_TOL)
   expect_true(any(grepl("favour the labelling in which component 1 is the minority", ws)))
  }
 }
 # weakly informative asymmetric priors: oriented, relabelled after sampling,
 # and no warning
 expect_no_warning(f_weak <- suppressMessages(glmMixBayes(X, y0, "gaussian",
                                                          priors = list(beta1 = "normal(0, 3)", beta2 = "normal(0, 4)"),
                                                          control = ctl)))
 expect_false(f_weak$diagnostics$orient_refused)
 expect_true(f_weak$diagnostics$relabelled)
 expect_gt(mean(f_weak$m_samples == 1L), 0.5)
})

test_that("predict.glmMixBayes validates newx and falls back to the stored design matrix", {
 d <- make_rev()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl, x = TRUE))
 X <- fit$x
 expect_length(predict(fit), nrow(d))
 expect_error(predict(fit, newx = X[1:3, 1, drop = FALSE]), "column")
 expect_error(predict(fit, newx = data.frame(X[1:3, ])), "numeric matrix")
 p0 <- predict(fit, newx = X[1:5, ])
 expect_equal(unname(p0), unname(as.vector(X[1:5, ] %*% colMeans(fit$estimates$coefficients))))
 p2 <- predict(fit, newx = X[1:5, ], interval = "credible", level = 0.9)
 expect_equal(colnames(p2), c("fit", "5 %", "95 %"))
 expect_true(all(p2[, 2] <= p2[, 3]))
 # without x = TRUE the design matrix is rebuilt from the stored model frame ...
 fit2 <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 expect_length(predict(fit2), nrow(d))
 # ... and a fit that stores neither needs `newx`
 fit3 <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl,
                                model = FALSE))
 expect_error(predict(fit3), "x = TRUE")

 expect_error(confint(fit, parm = "nope"), "Unknown parameter")
 expect_error(confint(fit, parm = 5), "indices")
 expect_equal(rownames(confint(fit, parm = 2)), "x")
 # parm = "theta" / "gamma" select the block when `block` is not given
 expect_equal(confint(fit, parm = "theta"), confint(fit, block = "theta"))
 expect_error(confint(fit, parm = "gamma"), "No 'gamma' draws")
})

test_that("survMixBayes methods cover the Weibull scale blocks and the Path B gamma block", {
 d <- make_rev()
 X <- cbind("(Intercept)" = 1, x = d$x)
 Z <- cbind("(Intercept)" = 1, z = d$z)
 S <- survival::Surv(d$time, d$status)
 fitW <- suppressMessages(survregMixBayes(X, S, dist = "weibull", control = ctl))
 s <- summary(fitW, probs = c(0.1, 0.9))
 expect_equal(names(confint(fitW)), c("coef1", "coef2", "theta", "shape1", "shape2", "scale1", "scale2",
                                      "m.coefficients"))
 expect_equal(colnames(s$coef1), c("Estimate", "Std. Error", "10 %", "90 %", "MCSE", "ESS", "Rhat"))
 expect_equal(rownames(s$coef1), c("(Intercept)", "x"))
 expect_output(print(s), "Scale multiplier of exp\\(X beta\\) \\(component 1 = correct-match")
 expect_output(print(s), "1 / shape (component 1 = correct-match; comparable to the Scale of survival::survreg())",
               fixed = TRUE)
 expect_equal(names(confint(fitW, parm = "scale1", level = 0.9)$scale1), c("5 %", "95 %"))
 pr <- predict(fitW, newdata = X[1:3, ], se.fit = TRUE, interval = "credible")
 expect_equal(colnames(pr$component1), c("fit", "se.fit", "2.5 %", "97.5 %"))
 expect_named(coef(fitW), c("(Intercept)", "x"))
 expect_error(predict(fitW, newdata = X[1:3, ], type = "response"), "Unused argument")
 expect_error(confint(fitW, parm = "nope"), "Unknown block")
 expect_error(confint(fitW, parm = "gamma"), "No 'gamma' draws")

 fitB <- suppressMessages(survregMixBayes(X, S, dist = "weibull", Z = Z, safe.matches = as.integer(d$safe),
                                          control = ctl))
 expect_output(print(summary(fitB)), "among records not flagged as safe matches")
 expect_equal(rownames(confint(fitB)$gamma), c("(Intercept)", "z"))
 expect_error(confint(fitB, parm = "theta"), "No 'theta' draws")

 fx <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "gamma",
                                  adjustment = adjMixBayes(d), control = ctl, x = TRUE))
 expect_length(predict(fx)$component1, nrow(d))
})

test_that("print.adjMixBayes() applies the prior precedence rules", {
 out <- capture_output(print(adjMixBayes(NULL, m.rate = 0.2)))
 expect_match(out, "theta ~ beta\\(12,3\\) \\[from m.rate\\]")
 expect_match(out, "Implied Prior on theta: median")
 out <- capture_output(print(adjMixBayes(NULL, m.rate = 0.2, priors = list(theta = "beta(2, 2)"))))
 expect_match(out, "not used: an explicit theta prior")
 expect_match(out, "theta: user-specified")
 df <- data.frame(z = stats::rnorm(10))
 out <- capture_output(print(adjMixBayes(df, m.formula = ~ z, priors = list(theta = "beta(2, 2)"))))
 expect_match(out, "logit moments of the theta prior")
 out <- capture_output(print(adjMixBayes(df, m.rate = 0.05)))
 expect_match(out, "density is unbounded")
 out <- capture_output(print(adjMixBayes(NULL, m.formula = ~ z, safe.matches = c(TRUE, FALSE))))
 expect_match(out, "Safe Matches:\\s+1 \\(50.0%\\)")
})
