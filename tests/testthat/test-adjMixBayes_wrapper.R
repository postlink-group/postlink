# End-to-end tests of the linkage information carried by adjMixBayes objects
# (m.formula, m.rate, safe.matches) through plglm() / plsurvreg() and into the
# C++ Gibbs sampler. Chains are short: these tests check wiring, output
# structure and gross posterior behaviour, not Monte Carlo precision.
local_edition(3)

make_linked <- function(n = 160, seed = 11) {
 set.seed(seed)
 x <- rnorm(n)
 z <- rnorm(n)                      # linkage covariate
 p_match <- stats::plogis(1.2 + 1.5 * z)
 match <- stats::rbinom(n, 1, p_match) == 1
 y <- ifelse(match, 1 + 2 * x + rnorm(n, sd = 0.7), rnorm(n, sd = 2))
 safe <- match & (seq_len(n) <= 40)   # the first 40 records: hand linked
 tt <- stats::rweibull(n, shape = 2, scale = exp(ifelse(match, 0.5 + 0.8 * x, 0.2)))
 cens <- stats::rexp(n, 0.1)
 data.frame(y = y, x = x, z = z, safe = safe,
            time = pmin(tt, cens), status = as.integer(tt <= cens))
}

ctrl <- list(iterations = 400, burnin.iterations = 200, seed = 5)

test_that("default adjMixBayes objects follow Path A (constant theta) through plglm()", {
 d <- make_linked()
 adj <- adjMixBayes(linked.data = d)
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adj, control = ctrl))

 expect_s3_class(fit, "glmMixBayes")
 expect_false(fit$use_logistic)
 expect_true(is.numeric(fit$estimates$theta))
 expect_length(fit$estimates$theta, 200)
 expect_null(fit$estimates$gamma)
 expect_equal(colnames(fit$estimates$coefficients), c("(Intercept)", "x"))

 s <- summary(fit)
 expect_false(s$use_logistic)
 expect_equal(rownames(s$theta), "theta")
 expect_equal(colnames(s$theta), c("Estimate", "Std. Error", "2.5 %", "97.5 %", "MCSE", "ESS", "Rhat"))
 expect_null(s$gamma)
 out <- capture_output(print(s))
 expect_match(out, "Match Probability theta (mixing weight of component 1 = correct-match):", fixed = TRUE)

 ci <- confint(fit, block = "theta")
 expect_equal(dim(ci), c(1L, 2L))
 expect_true(ci[1, 1] > 0 && ci[1, 2] < 1)
 expect_error(confint(fit, block = "gamma"), "No 'gamma' draws")
 expect_equal(rownames(confint(fit, block = "coefficients2")), c("(Intercept)", "x"))
 # the mismatch-indicator model: logit(1 - theta)
 expect_equal(rownames(confint(fit, block = "m.coefficients")), "(Intercept)")
 expect_equal(fit$estimates$m.coefficients[, "(Intercept)"], stats::qlogis(1 - fit$estimates$theta))
 expect_equal(rownames(confint(fit, parm = "x")), "x")
})

test_that("m.rate alone stays in Path A with a moment-matched beta prior", {
 d <- make_linked()
 adj <- adjMixBayes(linked.data = d, m.rate = 0.2, m.rate.sd = 0.05)
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adj, control = ctrl))

 expect_false(fit$use_logistic)
 expect_true(is.numeric(fit$estimates$theta))
 expect_null(fit$estimates$gamma)
 # informative prior centred at 0.8 pulls theta into a plausible range
 expect_true(abs(mean(fit$estimates$theta) - 0.8) < 0.15)
})

test_that("m.formula covariates activate Path B (gamma) and align Z with the outcome rows", {
 d <- make_linked()
 adj <- adjMixBayes(linked.data = d, m.formula = ~ z, m.rate = 0.2)
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adj, control = ctrl))

 expect_true(fit$use_logistic)
 expect_null(fit$estimates$theta)
 expect_true(is.matrix(fit$estimates$gamma))
 expect_equal(colnames(fit$estimates$gamma), c("(Intercept)", "z"))
 # z raises the match probability in the simulated data
 expect_gt(mean(fit$estimates$gamma[, "z"]), 0)

 s <- summary(fit)
 expect_true(s$use_logistic)
 expect_equal(rownames(s$gamma), c("(Intercept)", "z"))
 expect_null(s$theta)
 out <- capture_output(print(s))
 expect_match(out, "Match Probability Model gamma (logistic regression on the logit scale", fixed = TRUE)

 ci <- confint(fit, block = "gamma")
 expect_equal(rownames(ci), c("(Intercept)", "z"))
 expect_error(confint(fit, block = "theta"), "No 'theta' draws")

 # subsetting the outcome model keeps Z aligned through the row names
 d_na <- d
 d_na$y[c(3, 17, 90)] <- NA
 adj_na <- adjMixBayes(linked.data = d_na, m.formula = ~ z)
 fit_na <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adj_na, control = ctrl))
 expect_equal(ncol(fit_na$m_samples), nrow(d) - 3L)
 expect_equal(nrow(fit_na$estimates$gamma), 200)
})

test_that("safe.matches are fixed in component 1 and reported by print()", {
 d <- make_linked()
 adj <- adjMixBayes(linked.data = d, safe.matches = safe)
 expect_equal(sum(adj$safe.matches), sum(d$safe))

 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adj, control = ctrl))
 safe_cols <- which(d$safe)
 expect_true(all(fit$m_samples[, safe_cols] == 1L))
 expect_true(any(fit$m_samples[, -safe_cols] == 2L))

 out <- capture_output(print(adj))
 expect_match(out, sprintf("Safe Matches:\\s+%d \\(%.1f%%\\)", sum(d$safe), 100 * mean(d$safe)))
 expect_match(out, "Match Probability Model")

 # combined with covariates in m.formula
 adj2 <- adjMixBayes(linked.data = d, m.formula = ~ z, safe.matches = safe, m.rate = 0.2)
 fit2 <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adj2, control = ctrl))
 expect_true(all(fit2$m_samples[, safe_cols] == 1L))
 expect_equal(colnames(fit2$estimates$gamma), c("(Intercept)", "z"))
})

test_that("plsurvreg() carries m.formula and safe.matches into survregMixBayes()", {
 d <- make_linked()
 adj <- adjMixBayes(linked.data = d, m.formula = ~ z, safe.matches = safe)
 fit <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                  adjustment = adj, control = ctrl))

 expect_s3_class(fit, "survMixBayes")
 expect_true(fit$use_logistic)
 expect_equal(colnames(fit$estimates$gamma), c("(Intercept)", "z"))
 expect_equal(colnames(fit$estimates$coefficients), c("(Intercept)", "x"))
 expect_true(all(fit$m_samples[, which(d$safe)] == 1L))

 s <- summary(fit)
 expect_true(s$use_logistic)
 expect_equal(rownames(s$gamma), c("(Intercept)", "z"))
 expect_equal(rownames(s$m.coefficients), c("(Intercept)", "z"))
 expect_null(s$theta)
 out <- capture_output(print(s))
 expect_match(out, "Match Probability Model gamma (logistic regression on the logit scale", fixed = TRUE)
 expect_match(out, "gamma = -m.coefficients, among records not flagged as safe matches)", fixed = TRUE)

 ci <- confint(fit)
 expect_true("gamma" %in% names(ci))
 expect_false("theta" %in% names(ci))
 expect_equal(rownames(ci$gamma), c("(Intercept)", "z"))
 expect_equal(names(confint(fit, parm = "gamma")), "gamma")

 # Path A survival fits keep the theta block
 fitA <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                   adjustment = adjMixBayes(linked.data = d), control = ctrl))
 expect_false(fitA$use_logistic)
 expect_true("theta" %in% names(confint(fitA)))
 expect_true(!is.null(summary(fitA)$theta))
})

test_that("records with missing linkage covariates are dropped with a warning", {
 d <- make_linked()
 d$z[5] <- NA
 adj <- adjMixBayes(linked.data = d, m.formula = ~ z, safe.matches = safe)
 expect_warning(fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adj, control = ctrl)),
                "Dropped 1 observation")
 expect_equal(ncol(fit$m_samples), nrow(d) - 1L)
 expect_false("5" %in% fit$obs_names)
 expect_identical(colnames(fit$m_samples), fit$obs_names)
 expect_true(all(fit$m_samples[, as.character(which(d$safe))] == 1L))
 # infinite values from a transformation are an error
 d$z0 <- d$z; d$z0[is.na(d$z0)] <- 1; d$z0[7] <- 0
 expect_error(suppressWarnings(suppressMessages(plglm(y ~ x, family = "gaussian",
                                     adjustment = adjMixBayes(linked.data = d, m.formula = ~ log(z0)),
                                     control = ctrl))),
              "Infinite values")
})

test_that("the match probability model must keep its intercept", {
 d <- make_linked()
 expect_error(adjMixBayes(linked.data = d, m.formula = ~ z + 0), "must include an intercept")
 expect_error(adjMixBayes(linked.data = d, m.formula = ~ z - 1), "must include an intercept")
 # an adjustment object whose formula was altered afterwards is caught by the fitter
 adj <- adjMixBayes(linked.data = d, m.formula = ~ z)
 adj$m.formula <- ~ z - 1
 expect_error(suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adj, control = ctrl)),
              "must include an intercept")
})

test_that("user starting values are validated before reaching the sampler", {
 d <- make_linked()
 X <- cbind(1, d$x); colnames(X) <- c("(Intercept)", "x")
 small <- list(iterations = 60, burnin.iterations = 20, seed = 1)
 expect_error(suppressMessages(glmMixBayes(X, d$y, "gaussian", control = c(small, list(init = list(beta1 = 1))))),
              "length 2")
 expect_error(suppressMessages(glmMixBayes(X, d$y, "gaussian", control = c(small, list(init = list(foo = 1))))),
              "Unknown element")
 expect_error(suppressMessages(glmMixBayes(X, d$y, "gaussian", control = c(small, list(init = list(theta = NA))))),
              "finite numeric")
 fit <- suppressMessages(glmMixBayes(X, d$y, "gaussian",
                    control = c(small, list(init = list(beta1 = c(1, 2), beta2 = c(0, 0), theta = 0.7)))))
 expect_s3_class(fit, "glmMixBayes")
})
