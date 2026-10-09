# Names and methods of the Bayesian mixture fits aligned with the frequentist
# mixture fits of adjMixture(): the mismatch-indicator model in
# estimates$m.coefficients, component 2 in coefficients2 / dispersion2 /
# shape2 / scale2, the per-record match.prob, the summary columns and
# match-rate lines, the confint() and vcov() blocks, predict(newdata = ), the
# engine argument safe.matches and the generics that are refused.
local_edition(3)

make_align <- function(n = 160, seed = 11) {
 set.seed(seed)
 x <- stats::rnorm(n)
 z <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, stats::plogis(1.2 + 1.5 * z)) == 1
 y <- ifelse(match, 1 + 2 * x + stats::rnorm(n, sd = 0.7), stats::rnorm(n, sd = 2))
 tt <- stats::rweibull(n, shape = 2, scale = exp(ifelse(match, 0.5 + 0.8 * x, 0.2)))
 cens <- stats::rexp(n, 0.1)
 data.frame(y = y, x = x, f = factor(rep(c("a", "b", "c"), length.out = n)), z = z,
            safe = match & seq_len(n) <= 40,
            time = pmin(tt, cens), status = as.integer(tt <= cens))
}
ctl <- list(iterations = 400, burnin.iterations = 200, seed = 5)

test_that("estimates$m.coefficients is the mismatch-indicator model of glmMixture()", {
 d <- make_align()
 # Path A: a single column logit(1 - theta)
 fA <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 expect_equal(colnames(fA$estimates$m.coefficients), "(Intercept)")
 expect_equal(fA$estimates$m.coefficients[, "(Intercept)"], stats::qlogis(1 - fA$estimates$theta))
 expect_equal(dim(fA$estimates$coefficients2), dim(fA$estimates$coefficients))
 expect_equal(colnames(fA$estimates$coefficients2), c("(Intercept)", "x"))
 expect_true(all(c("dispersion", "dispersion2") %in% names(fA$estimates)))
 expect_false(any(c("m.dispersion", "m.shape", "m.scale") %in% names(fA$estimates)))

 # Path B: -gamma, named after the linkage covariates
 fB <- suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl,
                              adjustment = adjMixBayes(d, m.formula = ~ z, safe.matches = safe)))
 expect_equal(fB$estimates$m.coefficients, -fB$estimates$gamma)
 expect_equal(colnames(fB$estimates$m.coefficients), c("(Intercept)", "z"))
 # same sign as glmMixture(): z raises the match probability, so it lowers the mismatch probability
 expect_lt(mean(fB$estimates$m.coefficients[, "z"]), 0)

 # survival fits: the same layout with shape2 / scale2
 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull", control = ctl,
                                  adjustment = adjMixBayes(d, m.formula = ~ z)))
 expect_equal(fS$estimates$m.coefficients, -fS$estimates$gamma)
 expect_true(all(c("coefficients2", "shape2", "scale2") %in% names(fS$estimates)))
 fSA <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "gamma", control = ctl,
                                   adjustment = adjMixBayes(d)))
 expect_equal(fSA$estimates$m.coefficients[, 1], stats::qlogis(1 - fSA$estimates$theta))
 expect_null(fSA$estimates$scale2)
})

test_that("match.prob is the posterior share of draws in component 1, and 1 for safe matches", {
 d <- make_align()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl,
                               adjustment = adjMixBayes(d, safe.matches = safe)))
 expect_equal(fit$match.prob, colMeans(fit$m_samples == 1L))
 expect_identical(names(fit$match.prob), fit$obs_names)
 expect_true(all(fit$match.prob[d$safe] == 1))
 expect_true(all(fit$match.prob >= 0 & fit$match.prob <= 1))
 # placed after `call`, so the positions of the elements of 0.1.2 fits are unchanged
 expect_equal(names(fit)[1:5], c("m_samples", "estimates", "family", "call", "match.prob"))

 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull", control = ctl,
                                  adjustment = adjMixBayes(d, safe.matches = safe)))
 expect_equal(fS$match.prob, colMeans(fS$m_samples == 1L))
 expect_true(all(fS$match.prob[d$safe] == 1))

 # direct engine fits without row names give an unnamed vector
 X <- cbind(1, d$x)
 f0 <- suppressMessages(glmMixBayes(X, d$y, "gaussian", control = ctl))
 expect_null(names(f0$match.prob))
 expect_equal(unname(f0$match.prob), unname(colMeans(f0$m_samples == 1L)))
})

test_that("the engines take safe.matches (and the abbreviation safe)", {
 d <- make_align()
 X <- cbind("(Intercept)" = 1, x = d$x)
 f1 <- suppressMessages(glmMixBayes(X, d$y, "gaussian", safe.matches = d$safe, control = ctl))
 f2 <- suppressMessages(glmMixBayes(X, d$y, "gaussian", safe = d$safe, control = ctl))
 expect_identical(f1$m_samples, f2$m_samples)
 expect_equal(f1$diagnostics$n_safe, sum(d$safe))
 expect_error(suppressMessages(glmMixBayes(X, d$y, "gaussian", safe.matches = d$safe[-1], control = ctl)),
              "`safe.matches` must have length")
 # the sampler's own check uses the same name
 expect_error(postlink:::mixbayes_gibbs_cpp("gaussian", X, d$y, integer(0), matrix(0, nrow(X), 0),
                                            integer(1), list(), list(), 10L, 5L, 1L, 0L, 1L, TRUE, FALSE),
              "`safe.matches` must have length nrow\\(X\\)")
 S <- survival::Surv(d$time, d$status)
 s1 <- suppressMessages(survregMixBayes(X, S, "weibull", safe.matches = d$safe, control = ctl))
 s2 <- suppressMessages(survregMixBayes(X, S, "weibull", safe = d$safe, control = ctl))
 expect_identical(s1$m_samples, s2$m_samples)
 expect_true(all(s1$m_samples[, d$safe] == 1L))
 # through plglm() / plsurvreg() the linkage information comes from adjMixBayes()
 expect_error(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl,
                    safe.matches = d$safe), "given in adjMixBayes")
 expect_error(plsurvreg(survival::Surv(time, status) ~ x, adjustment = adjMixBayes(d), control = ctl,
                        m.rate = 0.1), "given in adjMixBayes")
})

test_that("summary() uses the confint() column names and prints both average match rates", {
 d <- make_align()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl,
                               adjustment = adjMixBayes(d, m.formula = ~ z, safe.matches = safe)))
 s <- summary(fit)
 cols <- c("Estimate", "Std. Error", "2.5 %", "97.5 %")
 # the tables of component 1 and of the match-probability model add the
 # Monte Carlo columns
 for (nm in c("coefficients", "m.coefficients", "gamma", "dispersion")) {
  expect_equal(colnames(s[[nm]]), c(cols, "MCSE", "ESS", "Rhat"), label = nm)
 }
 for (nm in c("coefficients2", "dispersion2")) expect_equal(colnames(s[[nm]]), cols, label = nm)
 expect_null(s$m.dispersion)
 expect_equal(rownames(s$m.coefficients), c("(Intercept)", "z"))
 expect_equal(s$m.coefficients[, "Estimate"], colMeans(fit$estimates$m.coefficients))
 expect_equal(s$m.coefficients[, "Std. Error"], apply(fit$estimates$m.coefficients, 2, stats::sd))
 expect_equal(s$coefficients2[, "Estimate"], colMeans(fit$estimates$coefficients2))

 # average correct-match probability over all records (safe ones count as 1)
 # and over the records not flagged as safe matches
 n <- nrow(d); n_safe <- sum(d$safe)
 expect_equal(s$match.rate$avg[["all"]], mean(fit$match.prob))
 expect_equal(s$match.rate$avg[["not.safe"]], mean(fit$match.prob[!d$safe]))
 expect_equal(s$match.rate$n, c(all = n, not.safe = n - n_safe))
 out <- paste(utils::capture.output(print(s)), collapse = "\n")
 # with safe matches every match-probability heading says, inside its
 # parentheses, that it describes the other records
 heads <- postlink:::.match_headings(TRUE)
 expect_match(out, paste0(heads[["m.coefficients"]], ":"), fixed = TRUE)
 expect_match(out, paste0(heads[["gamma"]], ":"), fixed = TRUE)
 expect_match(out, "Mismatch Model Coefficients (logit of the mismatch probability, among records not flagged as safe matches):",
              fixed = TRUE)
 expect_false(grepl("mismatch model describes", out))
 expect_match(out, sprintf("all records \\(n = %d; safe matches counted as 1\\): +%s", n,
                           format(signif(mean(fit$match.prob), 4))))
 expect_match(out, sprintf("records not flagged as safe matches \\(n = %d\\): +%s", n - n_safe,
                           format(signif(mean(fit$match.prob[!d$safe]), 4))))
 expect_match(out, "Std. Error: posterior standard deviation")

 # survival summary: released element names, the same columns
 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull", control = ctl,
                                  adjustment = adjMixBayes(d, safe.matches = safe)))
 sS <- summary(fS)
 for (nm in c("coef1", "m.coefficients", "theta", "shape1", "inv_shape1")) {
  expect_equal(colnames(sS[[nm]]), c(cols, "MCSE", "ESS", "Rhat"), label = nm)
 }
 for (nm in c("coef2", "shape2", "scale1", "scale2")) expect_equal(colnames(sS[[nm]]), cols, label = nm)
 expect_equal(rownames(sS$coef1), c("(Intercept)", "x"))
 expect_equal(unname(sS$coef1[, "Estimate"]), unname(colMeans(fS$estimates$coefficients)))
 expect_equal(colnames(summary(fS, probs = c(0.025, 0.5, 0.975))$coef1),
              c("Estimate", "Std. Error", "2.5 %", "50 %", "97.5 %", "MCSE", "ESS", "Rhat"))
 expect_error(summary(fS, probs = 2), "probabilities")
 expect_equal(sS$match.rate$avg[["not.safe"]], mean(fS$match.prob[!d$safe]))
 outS <- paste(utils::capture.output(print(sS)), collapse = "\n")
 # the same headings as the GLM summary
 expect_match(outS, paste0(heads[["m.coefficients"]], ":"), fixed = TRUE)
 expect_match(outS, paste0(heads[["theta"]], ":"), fixed = TRUE)
 expect_match(outS, "Match Probability theta (mixing weight of component 1 = correct-match, among records not flagged as safe matches):",
              fixed = TRUE)
 expect_false(grepl("Theta (mix weight", outS, fixed = TRUE))
 # without safe matches the headings carry no qualifier
 expect_equal(unname(postlink:::.match_headings(FALSE)[["theta"]]),
              "Match Probability theta (mixing weight of component 1 = correct-match)")
 expect_identical(postlink:::.match_headings(NULL), postlink:::.match_headings(FALSE))
 expect_match(outS, sprintf("records not flagged as safe matches \\(n = %d\\)", n - n_safe))
 expect_match(outS, "not survreg\\(\\)'s Scale")
})

test_that("confint() and vcov() select the same parameter blocks", {
 d <- make_align()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 est <- fit$estimates
 q95 <- function(v) unname(stats::quantile(v, c(0.025, 0.975)))
 expect_equal(colnames(confint(fit)), c("2.5 %", "97.5 %"))
 expect_equal(colnames(confint(fit, level = 0.9)), c("5 %", "95 %"))
 expect_equal(unname(confint(fit, block = "coefficients2")["x", ]), q95(est$coefficients2[, "x"]))
 expect_equal(unname(confint(fit, block = "m.coefficients")[1, ]), q95(est$m.coefficients[, 1]))
 expect_equal(unname(confint(fit, block = "dispersion2")[1, ]), q95(est$dispersion2))
 expect_equal(rownames(confint(fit, block = "dispersion")), "dispersion")
 # a single block name in `parm` selects the block, including "coefficients"
 expect_equal(confint(fit, parm = "coefficients2"), confint(fit, block = "coefficients2"))
 expect_equal(confint(fit, parm = "coefficients"), confint(fit))
 expect_equal(confint(fit, parm = "dispersion"), confint(fit, block = "dispersion"))
 expect_error(confint(fit, blok = "theta"), "Unused argument")

 for (b in c("coefficients", "coefficients2", "m.coefficients", "theta", "dispersion", "dispersion2")) {
  draws <- est[[b]]
  if (!is.matrix(draws)) draws <- matrix(draws, ncol = 1L, dimnames = list(NULL, b))
  expect_equal(vcov(fit, block = b), stats::cov(draws), label = b)
 }
 expect_equal(vcov(fit), stats::cov(est$coefficients))
 expect_error(vcov(fit, block = "gamma"), "No 'gamma' draws")
 expect_error(vcov(fit, blok = "theta"), "Unused argument")
 fP <- suppressMessages(glmMixBayes(cbind(1, rpois(160, 2)), rpois(160, 3), "poisson", control = ctl))
 expect_error(vcov(fP, block = "dispersion"), "no dispersion parameter")

 # survival: the list of 0.1.2 plus m.coefficients; parm with coefficient names or
 # indices and block = "<name>" give matrices
 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull", control = ctl,
                                  adjustment = adjMixBayes(d, m.formula = ~ z)))
 ci <- confint(fS)
 expect_true(is.list(ci))
 expect_equal(names(ci), c("coef1", "coef2", "gamma", "shape1", "shape2", "scale1", "scale2", "m.coefficients"))
 expect_equal(ci$m.coefficients, confint(fS, block = "m.coefficients"))
 cx <- confint(fS, parm = "x")
 expect_true(is.matrix(cx))
 expect_equal(cx, ci$coef1["x", , drop = FALSE])
 expect_equal(confint(fS, parm = 1:2), ci$coef1)
 expect_equal(confint(fS, block = "shape"), matrix(ci$shape1, 1L, dimnames = list("shape", names(ci$shape1))))
 expect_equal(confint(fS, block = "coefficients2", parm = "x"), ci$coef2["x", , drop = FALSE])
 expect_equal(names(confint(fS, parm = c("coef1", "gamma"))), c("coef1", "gamma"))
 expect_error(confint(fS, parm = c("x", "nope")), "Unknown block")
 # the block names accepted by `block` select the list elements in `parm` too
 expect_equal(confint(fS, parm = "coefficients"), ci["coef1"])
 expect_equal(confint(fS, parm = c("coefficients2", "shape", "scale", "m.coefficients")),
              ci[c("coef2", "shape1", "scale1", "m.coefficients")])
 expect_error(confint(fS, parm = c("shape", "nope")), "Unknown block\\(s\\) in `parm`: nope\\.")
 expect_error(confint(fS, parm = c("x", "gamma")), "separate calls")
 # unnamed columns (survregMixBayes() with an unnamed X) are named "(Intercept)"
 # and "X<j>"; coefficients can be selected by these names or by their indices
 fE <- suppressMessages(survregMixBayes(cbind(1, d$x), survival::Surv(d$time, d$status), "weibull",
                                        control = ctl))
 expect_equal(confint(fE, parm = "X2"), confint(fE, parm = 2))
 expect_equal(rownames(confint(fE, parm = 1:2)), c("(Intercept)", "X2"))
 expect_error(confint(fE, parm = "nope"), "indices of the coefficients of component 1: (Intercept), X2)",
              fixed = TRUE)
 expect_equal(vcov(fS, block = "m.coefficients"), stats::cov(fS$estimates$m.coefficients))
 expect_equal(vcov(fS, block = "scale2"), stats::cov(matrix(fS$estimates$scale2, ncol = 1L,
                                                            dimnames = list(NULL, "scale2"))))
 expect_error(vcov(fS, type = "x"), "Unused argument")
 # vcov() accepts the list names of confint() as its `block` argument does
 expect_equal(vcov(fS, block = "coef1"), vcov(fS))
 expect_equal(vcov(fS, block = "coef2"), vcov(fS, block = "coefficients2"))
 expect_equal(vcov(fS, block = "shape1"), vcov(fS, block = "shape"))
 expect_equal(vcov(fS, block = "scale1"), vcov(fS, block = "scale"))
 expect_error(vcov(fS, block = "nope"), "should be one of")
})

test_that("interval labels are formatted as by stats::confint() at every level", {
 d <- make_align()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 fl <- stats::lm(y ~ x, data = d)
 for (lv in c(0.5, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9999, 2 / 3)) {
  expect_equal(colnames(confint(fit, level = lv)), colnames(stats::confint(fl, level = lv)), label = lv)
 }
 expect_equal(colnames(confint(fit, level = 0.999)), c("0.05 %", "99.95 %"))
 expect_equal(colnames(predict(fit, newdata = d[1:2, ], interval = "credible", level = 0.995)),
              c("fit", "0.25 %", "99.75 %"))
 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                  adjustment = adjMixBayes(d), control = ctl))
 expect_equal(colnames(summary(fS, probs = c(0.0005, 0.5, 0.9995))$coef1),
              c("Estimate", "Std. Error", "0.05 %", "50 %", "99.95 %", "MCSE", "ESS", "Rhat"))
 expect_equal(colnames(confint(fS, level = 0.999)$coef1), c("0.05 %", "99.95 %"))
 expect_equal(colnames(predict(fS, newdata = cbind(1, d$x[1:2]), interval = "credible",
                               level = 0.999)$component1), c("fit", "0.05 %", "99.95 %"))
})

test_that("predict(newdata = ) builds the model matrix from the stored terms", {
 d <- make_align()
 fit <- suppressMessages(plglm(y ~ x + f, family = "gaussian", adjustment = adjMixBayes(d), control = ctl,
                               x = TRUE))
 nd <- d[c(5, 1, 9), ]
 expect_equal(predict(fit, newdata = nd), predict(fit, newx = fit$x[c(5, 1, 9), ]))
 expect_equal(predict(fit, newdata = nd, type = "response", interval = "credible"),
              predict(fit, newx = fit$x[c(5, 1, 9), ], type = "response", interval = "credible"))
 # factor levels come from the fit, even when the new data hold only one of them
 nd1 <- data.frame(x = c(0, 1), f = factor(c("b", "b")))
 X1 <- cbind("(Intercept)" = 1, x = c(0, 1), fb = 1, fc = 0)
 rownames(X1) <- rownames(nd1)
 expect_equal(predict(fit, newdata = nd1), predict(fit, newx = X1))
 # a numeric matrix in `newdata` is taken as the model matrix
 expect_equal(predict(fit, newdata = fit$x[1:3, ]), predict(fit, newx = fit$x[1:3, ]))
 expect_equal(colnames(predict(fit, newdata = nd, se.fit = TRUE, interval = "credible")),
              c("fit", "se.fit", "2.5 %", "97.5 %"))
 # missing covariates give NA rows (na.action = na.pass), also with intervals
 nd_na <- nd
 nd_na$x[2] <- NA
 p_na <- predict(fit, newdata = nd_na, se.fit = TRUE, interval = "credible")
 expect_true(all(is.na(p_na[2, ])))
 expect_false(anyNA(p_na[-2, ]))
 # misspelled or unused arguments are errors, not silently ignored
 expect_error(predict(fit, new_data = nd), "Unused argument")
 expect_error(predict(fit, nd, interval = "confidence"), "should be one of")
 expect_error(predict(fit, newdata = nd, newx = fit$x[1:3, ]), "either")
 expect_error(predict(fit, newdata = list(x = 1)), "data frame")

 # fits without stored terms need `newx`
 X <- cbind("(Intercept)" = 1, x = d$x)
 f0 <- suppressMessages(glmMixBayes(X, d$y, "gaussian", control = ctl))
 expect_error(predict(f0, newdata = d[1:3, ]), "pass the model matrix as `newx`")
 expect_length(predict(f0, newx = X[1:3, ]), 3L)

 # positional calls of postlink 0.1.2, predict(fit, X, type, se.fit, ...): the
 # second argument is newdata (a model matrix there), as in predict.glmMixture()
 expect_equal(names(formals(predict.glmMixBayes))[1:7], names(formals(predict.glmMixture))[1:7])
 expect_equal(predict(f0, X[1:3, ], "response"), predict(f0, newx = X[1:3, ], type = "response"))
 expect_equal(predict(f0, X[1:3, ], "link", TRUE), predict(f0, newx = X[1:3, ], se.fit = TRUE))
})

test_that("predict() without new data covers the analysed records only", {
 d <- make_align()
 d$z[c(3, 7)] <- NA   # records dropped at fit time (missing linkage covariate)
 fB <- suppressWarnings(suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl,
                                               adjustment = adjMixBayes(d, m.formula = ~ z))))
 fBx <- suppressWarnings(suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl, x = TRUE,
                                                adjustment = adjMixBayes(d, m.formula = ~ z))))
 n_used <- nrow(d) - 2L
 expect_equal(ncol(fB$m_samples), n_used)
 expect_length(predict(fB), n_used)          # from the stored model frame
 expect_length(predict(fBx), n_used)         # from the stored design matrix
 expect_equal(predict(fBx), predict(fBx, newx = fBx$x[fBx$obs_names, ]))
 expect_equal(predict(fB), predict(fB, newdata = d[-c(3, 7), ]))
 expect_equal(nrow(predict(fB, se.fit = TRUE)), n_used)

 fS <- suppressWarnings(suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                                   control = ctl,
                                                   adjustment = adjMixBayes(d, m.formula = ~ z))))
 expect_length(predict(fS)$component1, n_used)
 X_used <- cbind(1, d$x[-c(3, 7)])
 rownames(X_used) <- rownames(d)[-c(3, 7)]
 expect_equal(predict(fS), predict(fS, newdata = X_used))
})

test_that("predictions are named after the rows of the new data or the analysed records", {
 d <- make_align()
 d$z[c(3, 7)] <- NA
 fit <- suppressWarnings(suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl,
                                                adjustment = adjMixBayes(d, m.formula = ~ z))))
 nd <- d[c(9, 2), ]
 expect_equal(names(predict(fit, newdata = nd)), c("9", "2"))
 expect_equal(rownames(predict(fit, newdata = nd, type = "response", se.fit = TRUE, interval = "credible")),
              c("9", "2"))
 expect_equal(names(predict(fit)), names(fit$match.prob))      # the analysed records
 expect_equal(names(predict(fit)), rownames(d)[-c(3, 7)])
 expect_null(names(predict(fit, newx = cbind(1, c(0, 1)))))     # a matrix without row names
 fS <- suppressWarnings(suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                                   control = ctl,
                                                   adjustment = adjMixBayes(d, m.formula = ~ z))))
 expect_equal(names(predict(fS, newdata = nd)$component1), c("9", "2"))
 expect_equal(rownames(predict(fS, newdata = nd, se.fit = TRUE)$component2), c("9", "2"))
 expect_equal(names(predict(fS)$component2), names(fS$match.prob))
 expect_null(names(predict(fS, newdata = cbind(1, c(0, 1)))$component1))
})

test_that("predict() for survMixBayes fits takes a data frame in newdata", {
 d <- make_align()
 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x + f, dist = "weibull",
                                  adjustment = adjMixBayes(d), control = ctl))
 nd <- d[c(5, 1, 9), ]
 Xn <- stats::model.matrix(~ x + f, data = d)[c(5, 1, 9), ]
 expect_equal(predict(fS, newdata = nd), predict(fS, newdata = Xn))
 expect_equal(predict(fS, newdata = nd, se.fit = TRUE, interval = "credible"),
              predict(fS, newdata = Xn, se.fit = TRUE, interval = "credible"))
 # factor levels come from the fit
 nd1 <- data.frame(x = c(0, 1), f = factor(c("c", "c")))
 X1 <- cbind("(Intercept)" = 1, x = c(0, 1), fb = 0, fc = 1)
 rownames(X1) <- rownames(nd1)
 expect_equal(predict(fS, newdata = nd1), predict(fS, newdata = X1))
 # missing covariates give NA predictions (na.action = na.pass)
 nd2 <- nd
 nd2$x[2] <- NA
 expect_true(is.na(predict(fS, newdata = nd2)$component1[2]))
 p2 <- predict(fS, newdata = nd2, se.fit = TRUE, interval = "credible")
 for (k in c("component1", "component2")) {
  expect_true(all(is.na(p2[[k]][2, ])), label = k)
  expect_false(anyNA(p2[[k]][-2, ]), label = k)
 }
 # no new data: the records of the fit, from the stored model frame
 expect_equal(predict(fS), predict(fS, newdata = d))
 expect_error(predict(fS, newdata = nd, type = "response"), "Unused argument")
})

test_that("residual-based generics are refused for Bayesian fits only", {
 d <- make_align()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                  adjustment = adjMixBayes(d), control = ctl))
 for (obj in list(fit, fS)) {
  expect_error(stats::fitted(obj), "fitted\\(\\) is not available for Bayesian mixture fits")
  expect_error(stats::residuals(obj), "residuals\\(\\) is not available")
  expect_error(stats::resid(obj), "residuals\\(\\) is not available")
  expect_error(stats::df.residual(obj), "df.residual\\(\\) is not available")
  expect_error(stats::deviance(obj), "deviance\\(\\) is not available")
 }
 # the frequentist classes keep their methods
 mock <- structure(list(df.residual = 98, fitted.values = 1:3, residuals = c(0.1, -0.1, 0), deviance = 2),
                   class = c("glmMixture", "plglm", "plmodel"))
 expect_equal(stats::df.residual(mock), 98)
 expect_equal(stats::fitted(mock), 1:3)
 expect_equal(stats::residuals(mock), c(0.1, -0.1, 0))
 expect_equal(stats::deviance(mock), 2)
})

test_that("print() shows posterior means, the match model and the numbers of records", {
 d <- make_align()
 fA <- suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl,
                              adjustment = adjMixBayes(d, safe.matches = safe)))
 out <- paste(utils::capture.output(print(fA)), collapse = "\n")
 expect_match(out, "posterior means")
 expect_match(out, "Match probability theta")
 expect_match(out, "Dispersion \\(residual variance sigma\\^2")
 expect_match(out, sprintf("Records: %d \\(safe matches: %d\\)", nrow(d), sum(d$safe)))
 fB <- suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl,
                              adjustment = adjMixBayes(d, m.formula = ~ z)))
 expect_match(paste(utils::capture.output(print(fB)), collapse = "\n"), "Mismatch Model Coefficients")

 # survival: significant digits (small coefficients do not print as 0)
 X <- cbind("(Intercept)" = 1, x = d$x * 1e4)
 # (the covariate in large units may give a low acceptance rate in component 2)
 fS <- suppressWarnings(suppressMessages(survregMixBayes(X, survival::Surv(d$time, d$status), "weibull",
                                                         control = ctl)))
 outS <- paste(utils::capture.output(print(fS)), collapse = "\n")
 expect_match(outS, "e-0[45]")
 expect_match(outS, "Shape \\(component 1; posterior mean\\)")
 expect_match(outS, "Match probability theta")
})

test_that("mi_with() records and prints the refitted family and link", {
 d <- make_align()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 p <- mi_with(fit, data = d)
 expect_equal(c(p$refit, p$family, p$link), c("lm", "gaussian", "identity"))
 out <- paste(utils::capture.output(print(p)), collapse = "\n")
 expect_match(out, "Refit model: lm \\(family gaussian, link identity\\)")
 expect_match(out, "Estimate +Std. Error +2.5 % +97.5 % +df")
 # a Gaussian family with another link is refitted with glm()
 d$ypos <- exp(d$y / 4)
 pg <- mi_with(fit, data = d, formula = ypos ~ x, family = stats::gaussian(link = "log"))
 expect_equal(c(pg$refit, pg$family, pg$link), c("glm", "gaussian", "log"))
 expect_equal(mi_with(fit, data = d, family = "gaussian")$coef, p$coef)
 expect_error(mi_with(fit, data = d, family = 1), "family object")
 expect_error(mi_with(fit, data = d, family = "nofamily"), "Unknown family 'nofamily'")

 # the family name "gamma" (any case) is Gamma(link = "log"), as in glmMixture()
 # (base::gamma() is not a family); the family function Gamma keeps its inverse link
 d$g <- stats::rgamma(nrow(d), shape = 2, rate = 2 / exp(0.3 + 0.2 * d$x))
 fG <- suppressMessages(plglm(g ~ x, family = stats::Gamma(link = "log"), adjustment = adjMixBayes(d),
                              control = ctl))
 pG <- mi_with(fG, data = d)
 expect_equal(c(pG$refit, pG$family, pG$link), c("glm", "Gamma", "log"))
 for (nm in c("gamma", "Gamma", "GAMMA")) {
  pn <- mi_with(fG, data = d, family = nm)
  expect_equal(pn$link, "log", label = nm)
  expect_equal(pn$coef, pG$coef, label = nm)
 }
 expect_equal(mi_with(fG, data = d, family = stats::Gamma)$link, "inverse")

 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                  adjustment = adjMixBayes(d), control = ctl))
 outS <- paste(utils::capture.output(print(mi_with(fS, data = d))), collapse = "\n")
 expect_match(outS, "Estimate +Std. Error +2.5 % +97.5 % +df")
})

test_that("fits stored with the names of postlink 0.1.2 are read by the methods", {
 d <- make_align()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 old <- fit
 est <- fit$estimates
 old$estimates <- list(coefficients = est$coefficients, m.coefficients = est$coefficients2,
                       dispersion = est$dispersion, m.dispersion = est$dispersion2, theta = est$theta)
 old$match.prob <- NULL
 expect_equal(summary(old)$coefficients2, summary(fit)$coefficients2)
 expect_equal(summary(old)$m.coefficients, summary(fit)$m.coefficients)
 expect_equal(summary(old)$match.rate, summary(fit)$match.rate)
 expect_equal(confint(old, block = "m.coefficients"), confint(fit, block = "m.coefficients"))
 expect_equal(vcov(old, block = "dispersion2"), vcov(fit, block = "dispersion2"))
 expect_equal(posterior_draws(old), posterior_draws(fit))
 expect_output(print(old), "posterior means")

 # summary objects saved with the old names: the component-2 tables (then
 # m.coefficients, m.dispersion) are printed under component 2, and there is
 # no mismatch-model table to print
 s <- summary(fit)
 old_s <- s[c("call", "family", "coefficients", "dispersion", "use_logistic", "theta", "n_safe")]
 old_s$m.coefficients <- s$coefficients2
 old_s$m.dispersion <- s$dispersion2
 class(old_s) <- "summary.glmMixBayes"
 out_old <- utils::capture.output(res <- print(old_s))
 expect_identical(res, old_s)
 expect_false(any(grepl("NULL", out_old)))
 expect_false(any(grepl("Mismatch Model Coefficients", out_old)))
 comp2 <- grep("Component 2", out_old, fixed = TRUE)
 expect_length(comp2, 1L)
 expect_equal(out_old[comp2 + 1L], "Outcome Model Coefficients:")
 expect_match(out_old[comp2 + 3L], "^\\(Intercept\\)")
 expect_match(out_old[comp2 + 4L], "^x ")
 # the current summary prints the same component-2 table
 out_new <- utils::capture.output(print(s))
 k <- grep("Component 2", out_new, fixed = TRUE)
 expect_equal(out_old[comp2 + 0:4], out_new[k + 0:4])

 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                  adjustment = adjMixBayes(d), control = ctl))
 oldS <- fS
 e <- fS$estimates
 oldS$estimates <- list(coefficients = e$coefficients, m.coefficients = e$coefficients2, theta = e$theta,
                        shape = e$shape, m.shape = e$shape2, scale = e$scale, m.scale = e$scale2)
 expect_equal(predict(oldS, newdata = cbind(1, d$x[1:3])), predict(fS, newdata = cbind(1, d$x[1:3])))
 expect_equal(confint(oldS), confint(fS))
 expect_equal(summary(oldS)$scale2, summary(fS)$scale2)
})
