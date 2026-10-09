# Warnings and diagnostics of the Bayesian mixture fits: the prior-scale and
# prior-data conflict check, the effective sample size and split R-hat checks,
# the check of the log posterior, the window of the acceptance rates, the
# Monte Carlo columns and notes of summary(), the default burn-in and the
# minimum number of draws of mi_with().
local_edition(3)

make_warn <- function(n = 150, seed = 3) {
 set.seed(seed)
 x <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, 0.8) == 1
 y <- ifelse(match, 1 + 2 * x + stats::rnorm(n, sd = 0.7), stats::rnorm(n, sd = 2))
 tt <- stats::rweibull(n, shape = 2, scale = exp(ifelse(match, 0.5 + 0.8 * x, 0.2)))
 data.frame(y = y, x = x, z = stats::rnorm(n), time = tt, status = 1L)
}
ctl <- list(iterations = 400, burnin.iterations = 200, seed = 4)

# warnings of an expression, which is evaluated with messages suppressed
collect_warnings <- function(expr) {
 ws <- character()
 value <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
  ws <<- c(ws, conditionMessage(w))
  invokeRestart("muffleWarning")
 })
 list(value = value, warnings = ws)
}

# ------------------------------------------------------------------------------
# Prior-data conflict and prior scale
# ------------------------------------------------------------------------------
test_that("the prior-data conflict warning says which priors are defaults and names the help page", {
 b1 <- cbind(rep(50, 10), rep(0.1, 10))
 ep <- list(beta1_mean = c(0, 0), beta1_sd = c(10, 5), intercept = TRUE)
 expect_warning(postlink:::.check_prior_data_conflict(b1, ep),
                "default prior of '\\(Intercept\\)'.*not rescaled to the data: they suit outcomes and covariates of order one")
 w <- tryCatch(postlink:::.check_prior_data_conflict(b1, ep, default = c(FALSE, TRUE)),
               warning = conditionMessage)
 expect_match(w, "supplied prior of '(Intercept)'", fixed = TRUE)
 expect_match(w, "Check that the supplied priors suit the scale", fixed = TRUE)
 expect_false(grepl("default", w))
 # survival fits refer to ?survregMixBayes
 expect_warning(postlink:::.check_prior_data_conflict(b1, ep, help = "survregMixBayes"),
                "(see 'Prior distributions' in ?survregMixBayes).", fixed = TRUE)
 # unnamed design columns: an intercept column is "(Intercept)", others "X<j>"
 b2 <- cbind(rep(0, 10), rep(50, 10))
 expect_warning(postlink:::.check_prior_data_conflict(b2, ep), "default prior of 'X2'")
 expect_equal(postlink:::.coef_labels(c("", "x", NA), 3L, intercept = TRUE), c("(Intercept)", "x", "X3"))
 expect_equal(postlink:::.coef_labels(NULL, 2L, intercept = FALSE), c("X1", "X2"))
 expect_false(postlink:::.check_prior_data_conflict(cbind(rep(1, 10), 0), ep))

 # the naive fit is described by its source: only a fit to all records
 # ignores the linkage errors
 set.seed(1)
 b0 <- cbind(stats::rnorm(200, 0.5, 10), stats::rnorm(200, 0, 5))
 nv_safe <- list(est = c(30000, 2000), source = "the safe matches")
 w <- tryCatch(postlink:::.check_prior_data_conflict(b0, ep, naive = nv_safe), warning = conditionMessage)
 expect_match(w, "a fit of the outcome model to the safe matches gives 30000, 3000.0 prior standard deviations",
              fixed = TRUE)
 expect_false(grepl("ignores the linkage errors", w))
 expect_match(w, "stays close to the prior: the prior dominates the estimate", fixed = TRUE)
 w <- tryCatch(postlink:::.check_prior_data_conflict(b0, ep, naive = list(est = c(30000, 2000), source = "all records")),
               warning = conditionMessage)
 expect_match(w, "a fit of the outcome model to all records that ignores the linkage errors gives 30000", fixed = TRUE)

 # a posterior whose standard deviation stays close to the prior standard
 # deviation but whose mean lies far from the prior mean: the prior still
 # dominates the estimate, and the message does not say that the posterior
 # stays close to the prior
 b_far <- cbind(stats::rnorm(200, 0, 10), stats::rnorm(200, 22, 5))
 w <- tryCatch(postlink:::.check_prior_data_conflict(b_far, ep, naive = list(est = c(0, 300), source = "all records")),
               warning = conditionMessage)
 expect_match(w, "The default prior of 'X2'", fixed = TRUE)
 expect_match(w, paste0("and the posterior mean \\([0-9.]+\\) lies 4\\.[0-9] prior standard deviations from it, while ",
                        "the posterior standard deviation \\([0-9.]+\\) stays close to the prior standard deviation: ",
                        "the prior still dominates the estimate"))
 expect_false(grepl("stays close to the prior:", w, fixed = TRUE))

 # a supplied prior is flagged after sampling when it dominates the estimate
 # (posterior SD above 0.9 prior SDs) although the naive estimate contradicts
 # it, also allowing for its standard error
 set.seed(2)
 ep2 <- list(beta1_mean = c(0, 0), beta1_sd = c(10, 2), intercept = TRUE)
 b_dom <- cbind(stats::rnorm(200, 0, 10), 0.3 + 1.95 * as.numeric(scale(stats::rnorm(200))))
 nv_far <- list(est = c(0.5, 495), se = c(0.1, 40), source = "all records")
 w <- tryCatch(postlink:::.check_prior_data_conflict(b_dom, ep2, default = c(TRUE, FALSE), naive = nv_far),
               warning = conditionMessage)
 expect_match(w, paste0("The supplied prior of 'X2' in the correct-match component, normal(0, 2), conflicts with ",
                        "the data: a fit of the outcome model to all records that ignores the linkage ",
                        "errors gives 495, 247.5 prior standard deviations from the prior mean, while the posterior ",
                        "(mean 0.3, standard deviation 1.95) stays close to the prior: the prior dominates the estimate."),
              fixed = TRUE)
 expect_match(w, "Check that the supplied priors suit the scale of the outcome and covariates and are not more informative than intended",
              fixed = TRUE)
 # weak data (a naive estimate within 3 standard deviations of the prior mean
 # once its standard error is allowed for) may be dominated by a supplied
 # informative prior without a warning, but not by a default prior
 nv_weak <- list(est = c(0.5, 7), se = c(0.1, 3), source = "all records")
 expect_false(postlink:::.check_prior_data_conflict(b_dom, ep2, default = c(TRUE, FALSE), naive = nv_weak))
 expect_warning(postlink:::.check_prior_data_conflict(b_dom, ep2, naive = nv_weak), "The default prior of 'X2'")
 # without a standard error, or when the prior does not dominate, or before
 # sampling, a supplied prior is not flagged by the naive fit
 expect_false(postlink:::.check_prior_data_conflict(b_dom, ep2, default = c(TRUE, FALSE),
                                                    naive = list(est = c(0.5, 495), source = "all records")))
 b_ok <- cbind(stats::rnorm(200, 0, 10), stats::rnorm(200, 0.3, 0.5))
 expect_false(postlink:::.check_prior_data_conflict(b_ok, ep2, default = c(TRUE, FALSE), naive = nv_far))
 expect_false(postlink:::.check_prior_data_conflict(NULL, ep2, default = c(TRUE, FALSE), naive = nv_far,
                                                    coef_names = c("(Intercept)", "x")))

 # Weibull fits: the posterior mean of the identified intercept (with
 # log(scale) added, as summary() reports it) is quoted next to the sampled one
 b_w <- cbind(43 + 9.5 * as.numeric(scale(stats::rnorm(200))), stats::rnorm(200, 0, 0.1))
 w <- tryCatch(postlink:::.check_prior_data_conflict(b_w, ep2, help = "survregMixBayes",
                                                     naive = list(est = c(691, 0), se = c(1, 0.1), source = "all records"),
                                                     coef_names = c("(Intercept)", "x"),
                                                     reported = c(43.6, 0)),
               warning = conditionMessage)
 expect_match(w, paste0("the posterior mean (43) lies 4.3 prior standard deviations from it, while the posterior ",
                        "standard deviation (9.5) stays close to the prior standard deviation: the prior still ",
                        "dominates the estimate (the prior and the posterior quoted here refer to the sampled ",
                        "intercept; summary() reports the intercept with log(scale) added, whose posterior mean ",
                        "is 43.6)."), fixed = TRUE)

 # when sampling stopped (no draws), only the scale check is reported, with
 # huge distances in exponent notation
 w <- tryCatch(postlink:::.check_prior_data_conflict(NULL, ep, naive = list(est = c(1e201, 1e200), source = "all records"),
                                                     coef_names = c("", "x")),
               warning = conditionMessage)
 expect_match(w, "The default prior of '(Intercept)' in the correct-match component, normal(0, 10), does not suit the scale of the data",
              fixed = TRUE)
 expect_match(w, "gives 1e+201, 1e+200 prior standard deviations from the prior mean (2 coefficients are affected; also 'x').",
              fixed = TRUE)
 expect_false(grepl("posterior", w))

 # which coefficients use a default prior
 expect_equal(postlink:::.default_prior_columns(NULL, 3L, TRUE), rep(TRUE, 3))
 expect_equal(postlink:::.default_prior_columns(list(intercept1 = "normal(0, 1)"), 3L, TRUE), c(FALSE, TRUE, TRUE))
 expect_equal(postlink:::.default_prior_columns(list(beta1 = "normal(0, 1)"), 3L, TRUE), c(TRUE, FALSE, FALSE))
 expect_equal(postlink:::.default_prior_columns(list(intercept1 = "normal(0, 1)"), 2L, FALSE), c(TRUE, TRUE))
 # an entry that is present but NULL is not supplied (as in prepare_mixbayes_priors())
 expect_equal(postlink:::.default_prior_columns(list(beta1 = NULL, intercept1 = "normal(0, 1)"), 3L, TRUE),
              c(FALSE, TRUE, TRUE))
 expect_equal(postlink:::.default_prior_columns(list(beta1 = NULL, intercept1 = NULL), 2L, TRUE), c(TRUE, TRUE))

 # end to end: a prior supplied by the user is not called a default, even
 # when it equals the default (the outcome is on the scale of years)
 d <- make_warn()
 d$age <- 60 + 5 * d$y
 r <- collect_warnings(plglm(age ~ x, family = "gaussian", control = ctl,
                             adjustment = adjMixBayes(d, priors = list(intercept1 = "normal(0, 10)",
                                                                       beta1 = "normal(0, 20)"))))
 pw <- grep("prior of", r$warnings, value = TRUE)
 expect_length(pw, 1L)
 expect_match(pw, "The supplied prior of '(Intercept)' in the correct-match component, normal(0, 10), is far from the data",
              fixed = TRUE)
 expect_false(grepl("default", pw))
 # with the default intercept prior the same fit names it a default
 r0 <- collect_warnings(plglm(age ~ x, family = "gaussian", control = ctl,
                              adjustment = adjMixBayes(d, priors = list(beta1 = "normal(0, 20)"))))
 expect_match(grep("prior of", r0$warnings, value = TRUE), "The default prior of '(Intercept)'", fixed = TRUE)
})

test_that("default priors that do not suit the scale of the data are reported, once per fit", {
 # a Gaussian outcome in very large units (CRIT-1): the residual SD absorbs
 # the outcome and the posterior collapses onto the default priors, so the
 # posterior means stay within 3 prior standard deviations
 set.seed(10)
 n <- 300
 x <- stats::rnorm(n)
 ok <- stats::runif(n) < 0.9
 base <- 3.3 + 0.2 * x + 0.5 * stats::rnorm(n)
 base[!ok] <- sample(base[ok], sum(!ok), replace = TRUE)
 X <- cbind("(Intercept)" = 1, x = x)
 ctl_s <- list(iterations = 600, burnin.iterations = 300, seed = 1)
 r <- collect_warnings(glmMixBayes(X, 1e4 * base, "gaussian", control = ctl_s))
 pw <- grep("prior of", r$warnings, value = TRUE)
 expect_length(pw, 1L)
 expect_match(pw, "The default prior of '(Intercept)' in the correct-match component, normal(0, 10), does not suit the scale of the data",
              fixed = TRUE)
 expect_match(pw, "ignores the linkage errors gives 3[0-9]{4}, [0-9.]+ prior standard deviations from the prior mean")
 expect_match(pw, "the prior dominates the estimate (2 coefficients are affected; also 'x')", fixed = TRUE)
 # the posterior means alone are not flagged
 ep <- list(beta1_mean = c(0, 0), beta1_sd = c(10, 5), intercept = TRUE)
 expect_false(postlink:::.check_prior_data_conflict(r$value$estimates$coefficients, ep))

 # the same priors supplied by the user are flagged as supplied priors that
 # dominate the estimate
 r1 <- collect_warnings(glmMixBayes(X, 1e4 * base, "gaussian", control = ctl_s,
                                    priors = list(intercept1 = "normal(0, 10)", beta1 = "normal(0, 5)")))
 pw1 <- grep("prior of", r1$warnings, value = TRUE)
 expect_length(pw1, 1L)
 expect_match(pw1, "The supplied prior of '(Intercept)' in the correct-match component, normal(0, 10), conflicts with the data",
              fixed = TRUE)
 expect_match(pw1, "the prior dominates the estimate (2 coefficients are affected; also 'x')", fixed = TRUE)
 expect_false(grepl("default", pw1))

 # priors on the scale of the data: no warning on the priors
 scaled <- list(intercept1 = "normal(30000, 10000)", beta1 = "normal(0, 10000)",
                intercept2 = "normal(30000, 10000)", beta2 = "normal(0, 10000)",
                sigma1 = "cauchy(0, 5000)", sigma2 = "cauchy(0, 5000)")
 r2 <- collect_warnings(glmMixBayes(X, 1e4 * base, "gaussian", priors = scaled, control = ctl_s))
 expect_false(any(grepl("prior of", r2$warnings)))

 # a covariate in very small units in a Weibull model: the survival fit refers
 # to ?survregMixBayes
 tt <- stats::rweibull(n, 1.5, exp(0.5 + 0.8 * x))
 tt[!ok] <- sample(tt[ok], sum(!ok), TRUE)
 Xs <- cbind("(Intercept)" = 1, xs = x / 1000)
 r3 <- collect_warnings(survregMixBayes(Xs, cbind(tt, 1), "weibull", control = ctl_s))
 pw3 <- grep("prior of", r3$warnings, value = TRUE)
 expect_length(pw3, 1L)
 expect_match(pw3, "The default prior of 'xs' in the correct-match component, normal(0, 2), does not suit", fixed = TRUE)
 expect_match(pw3, "?survregMixBayes", fixed = TRUE)
 # the "stays close to the prior" wording requires a posterior mean within 3
 # prior standard deviations of the prior mean
 b <- r3$value$estimates$coefficients[, "xs"]
 if (abs(mean(b)) > 3 * 2) {
  expect_false(grepl("stays close to the prior:", pw3, fixed = TRUE))
  expect_match(pw3, "the prior still dominates the estimate|lies [0-9.]+ prior standard deviations from it")
 } else {
  expect_false(grepl("still dominates", pw3, fixed = TRUE))
 }
 # the flagged coefficient is not the intercept: no note on log(scale)
 expect_false(grepl("log(scale)", pw3, fixed = TRUE))

 # survival times in very large units: the Weibull intercept is flagged, and
 # the message quotes the posterior mean of the intercept that summary()
 # reports (with log(scale) added)
 d <- make_warn()
 Xw <- cbind("(Intercept)" = 1, x = d$x)
 r4 <- collect_warnings(survregMixBayes(Xw, cbind(d$time * 1e30, 1), "weibull", control = ctl))
 pw4 <- grep("prior of", r4$warnings, value = TRUE)
 expect_length(pw4, 1L)
 expect_match(pw4, "The default prior of '(Intercept)'", fixed = TRUE)
 expect_match(pw4, sprintf("summary() reports the intercept with log(scale) added, whose posterior mean is %s)",
                           format(signif(summary(r4$value)$coef1["(Intercept)", "Estimate"], 4))), fixed = TRUE)
 # the gamma distribution has no scale term
 r5 <- collect_warnings(survregMixBayes(Xw, cbind(d$time * 1e30, 1), "gamma", control = ctl))
 expect_false(any(grepl("log(scale)", r5$warnings, fixed = TRUE)))
})

test_that("the naive fit of the prior-scale check uses the safe matches and gives up on warnings", {
 set.seed(5)
 n <- 100
 x <- stats::rnorm(n)
 X <- cbind("(Intercept)" = 1, x = x)
 y <- 1 + 2 * x + stats::rnorm(n)
 nv <- postlink:::.naive_coefficients(X, y, "gaussian")
 expect_equal(nv$source, "all records")
 expect_equal(nv$est, unname(stats::coef(stats::lm(y ~ x))))
 expect_equal(nv$se, unname(summary(stats::lm(y ~ x))$coefficients[, "Std. Error"]))
 # with at least p + 5 safe matches only these are used
 safe <- c(rep(1L, 7), rep(0L, n - 7))
 nv7 <- postlink:::.naive_coefficients(X, y, "gaussian", safe = safe)
 expect_equal(nv7$source, "the safe matches")
 expect_equal(nv7$est, unname(stats::coef(stats::lm(y ~ x, subset = safe == 1L))))
 expect_equal(postlink:::.naive_coefficients(X, y, "gaussian", safe = c(rep(1L, 6), rep(0L, n - 6)))$source,
              "all records")
 # GLM families use the links of the mixture components
 yp <- stats::rpois(n, exp(0.5 + 0.3 * x))
 expect_equal(postlink:::.naive_coefficients(X, yp, "poisson")$est,
              unname(stats::coef(stats::glm(yp ~ x, family = stats::poisson()))))
 expect_equal(postlink:::.naive_coefficients(X, yp, "poisson")$se,
              unname(summary(stats::glm(yp ~ x, family = stats::poisson()))$coefficients[, "Std. Error"]))
 yg <- stats::rgamma(n, 2, rate = 2 / exp(0.5 + 0.3 * x))
 nvg <- postlink:::.naive_coefficients(X, yg, "gamma")
 gg <- stats::glm(yg ~ x, family = stats::Gamma(link = "log"))
 expect_equal(nvg$est, unname(stats::coef(gg)), tolerance = 1e-6)
 expect_equal(nvg$se, unname(summary(gg)$coefficients[, "Std. Error"]), tolerance = 1e-6)
 # a fit that warns (here: complete separation) is not used
 yb <- as.integer(x > 0)
 expect_null(postlink:::.naive_coefficients(X, yb, "binomial"))
 # survival models: a Weibull accelerated failure time fit
 tt <- stats::rweibull(n, 1.5, exp(0.5 + 0.8 * x))
 ev <- rep(1L, n)
 sr <- survival::survreg(survival::Surv(tt, ev) ~ x, dist = "weibull")
 nvs <- postlink:::.naive_coefficients(X, tt, "weibull", event = ev)
 expect_equal(nvs$est, unname(stats::coef(sr)), tolerance = 1e-6)
 expect_equal(nvs$se, unname(sqrt(diag(stats::vcov(sr)))[1:2]), tolerance = 1e-5)
 expect_null(postlink:::.naive_coefficients(X, tt, "weibull", event = rep(0L, n)))

 # when the fit to the safe matches fails or warns, all records are used:
 # no events among the safe matches, an outcome without variation among them
 ev_safe0 <- ifelse(safe == 1L, 0L, 1L)
 nvf <- postlink:::.naive_coefficients(X, tt, "weibull", event = ev_safe0, safe = safe)
 expect_equal(nvf$source, "all records")
 expect_equal(nvf$est, unname(stats::coef(survival::survreg(survival::Surv(tt, ev_safe0) ~ x, dist = "weibull"))),
              tolerance = 1e-6)
 yb2 <- stats::rbinom(n, 1, stats::plogis(x))
 yb2[safe == 1L] <- 1L
 nvb <- postlink:::.naive_coefficients(X, yb2, "binomial", safe = safe)
 expect_equal(nvb$source, "all records")
 expect_equal(nvb$est, unname(stats::coef(stats::glm(yb2 ~ x, family = stats::binomial()))))
 # and nothing is concluded when that fit fails too
 expect_null(postlink:::.naive_coefficients(X, rep(1L, n), "binomial", safe = safe))
})

# ------------------------------------------------------------------------------
# Effective sample size, split R-hat and the log posterior
# ------------------------------------------------------------------------------
test_that("the low-ESS warning is gated on the iterations after the burn-in, not on the stored draws", {
 low <- c("coefficients[x]" = 40, theta = 300)
 expect_warning(postlink:::.check_ess(low, n_draws = 250, n_sweeps = 5000),
                "coefficients[x] 40 (out of 250 stored draws from 5000 iterations after the burn-in)",
                fixed = TRUE)
 expect_silent(postlink:::.check_ess(low, n_draws = 3000, n_sweeps = 999))
 # a thinned run storing fewer than 200 draws is also advised to reduce `thin`
 expect_warning(postlink:::.check_ess(low, n_draws = 100, n_sweeps = 10000), "or reduce `thin`", fixed = TRUE)
 w <- tryCatch(postlink:::.check_ess(low, n_draws = 250, n_sweeps = 5000), warning = conditionMessage)
 expect_false(grepl("thin", w))
 w <- tryCatch(postlink:::.check_ess(low, n_draws = 1000, n_sweeps = 1000), warning = conditionMessage)
 expect_false(grepl("thin", w))

 # thinning does not silence the warning
 local_mocked_bindings(.mixbayes_ess = function(fit) {
  d <- postlink:::.mixbayes_draws(fit)
  stats::setNames(ifelse(colnames(d) == "theta", 20, 500), colnames(d))
 }, .package = "postlink")
 d <- make_warn()
 X <- cbind("(Intercept)" = 1, x = d$x)
 r <- collect_warnings(glmMixBayes(X, d$y, "gaussian",
                                    control = list(iterations = 1300, burnin.iterations = 300, thin = 10,
                                                   seed = 1)))
 expect_true(any(grepl("theta 20 (out of 100 stored draws from 1000 iterations after the burn-in)",
                       r$warnings, fixed = TRUE)))
 r2 <- collect_warnings(glmMixBayes(X, d$y, "gaussian",
                                    control = list(iterations = 1290, burnin.iterations = 300, seed = 1)))
 expect_false(any(grepl("Low effective sample size", r2$warnings)))
})

test_that("the split R-hat of one chain equals rstan::Rhat()", {
 # reference values computed with rstan 2.32.7, rstan::Rhat(x) (equal to
 # rstan::Rhat(matrix(x, ncol = 1)), the draws reshaped as one chain); they
 # are hard-coded so that the tests do not depend on rstan
 set.seed(3)
 xs <- list(stats::rnorm(1000), cumsum(stats::rnorm(1001)), stats::rexp(57)^2,
            c(stats::rnorm(500), stats::rnorm(500, 1)), stats::rt(2000, 3))
 ref <- c(1.00114234546701, 1.12829224199753, 1.01652950968873, 1.22012099989854, 1.00233356298913)
 expect_equal(vapply(xs, postlink:::.rhat_split, numeric(1)), ref, tolerance = 1e-12)
 # ties, and draws on a tiny scale: constancy is judged after the rank
 # transformation, as in rstan
 set.seed(1)
 tiny <- 1e-18 * stats::rnorm(200)
 ties <- round(cumsum(stats::rnorm(150)), 0)
 expect_equal(postlink:::.rhat_split(tiny), 1.00185886745232, tolerance = 1e-12)
 expect_equal(postlink:::.rhat_split(ties), 1.4097869405414, tolerance = 1e-12)
 # a chain stuck at a different value in each half: Inf (rstan gives NA,
 # since the folded draws are constant), and it is warned about
 step <- c(rep(1, 50), rep(2, 50))
 expect_equal(postlink:::.rhat_split(step), Inf)
 expect_equal(postlink:::.high_rhat(c(theta = Inf, lp = NA, x = 1.01)), c(theta = Inf))
 expect_warning(postlink:::.check_rhat(c(theta = Inf), n_draws = 1000, n_sweeps = 1000), "theta Inf")
})

test_that("split R-hat is stored for the reported parameters, warned about and noted by summary()", {
 # basic properties: NA for constant or too short chains, large for drift
 expect_true(is.na(postlink:::.rhat_split(rep(2, 100))))
 expect_true(is.na(postlink:::.rhat_split(1:3)))
 expect_true(is.na(postlink:::.rhat_split(c(1, 2, NA, 4, 5))))
 set.seed(1)
 expect_lt(postlink:::.rhat_split(stats::rnorm(2000)), 1.01)
 expect_gt(postlink:::.rhat_split(c(stats::rnorm(500), stats::rnorm(500, 3))), 1.5)

 d <- make_warn()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 rh <- fit$diagnostics$rhat
 expect_named(rh, c("coefficients[(Intercept)]", "coefficients[x]", "theta", "dispersion", "lp"))
 draws <- posterior_draws(fit)
 expect_equal(unname(rh), unname(apply(draws[, names(rh)], 2, postlink:::.rhat_split)))
 # with linkage covariates: the mismatch-model coefficients; survival: the shape
 fB <- suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl,
                              adjustment = adjMixBayes(d, m.formula = ~ z)))
 expect_named(fB$diagnostics$rhat, c("coefficients[(Intercept)]", "coefficients[x]",
                                     "m.coefficients[(Intercept)]", "m.coefficients[z]", "dispersion", "lp"))
 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                  adjustment = adjMixBayes(d), control = ctl))
 expect_named(fS$diagnostics$rhat, c("coefficients[(Intercept)]", "coefficients[x]", "theta", "shape", "lp"))

 # the warning applies to runs with at least 1000 iterations after the burn-in
 # and at least 200 stored draws (heavily thinned runs give false alarms)
 expect_warning(postlink:::.check_rhat(c(theta = 1.2, lp = 1.01), n_draws = 500, n_sweeps = 1000),
                "Split R-hat above 1.05: theta 1.200 (single chain of 500 stored draws)", fixed = TRUE)
 expect_warning(postlink:::.check_rhat(c(theta = 1.2), n_draws = 500, n_sweeps = 1000), "several seeds")
 expect_silent(postlink:::.check_rhat(c(theta = 1.2), n_draws = 500, n_sweeps = 999))
 expect_silent(postlink:::.check_rhat(c(theta = 1.01, lp = NA), n_draws = 5000, n_sweeps = 5000))
 expect_warning(postlink:::.check_rhat(c(theta = 1.2), n_draws = 200, n_sweeps = 20000), "theta 1.200")
 expect_silent(postlink:::.check_rhat(c(theta = 1.2), n_draws = 199, n_sweeps = 20000))

 # summary() notes high values for any run
 fit$diagnostics$rhat[] <- 1
 expect_length(summary(fit)$high_rhat, 0L)
 fit$diagnostics$rhat[["theta"]] <- 1.3
 s <- summary(fit)
 expect_equal(s$high_rhat, c(theta = 1.3))
 expect_equal(unname(s$theta[, "Rhat"]), 1.3)
 expect_match(paste(utils::capture.output(print(s)), collapse = "\n"),
              "Note: split R-hat above 1.05 for theta (1.300)", fixed = TRUE)
 fS$diagnostics$rhat[] <- 1
 fS$diagnostics$rhat[["shape"]] <- 1.2
 expect_match(paste(utils::capture.output(print(summary(fS))), collapse = "\n"),
              "split R-hat above 1.05 for shape (1.200);", fixed = TRUE)
 # the note says that R-hat is noisy when fewer than 200 draws were stored
 expect_equal(s$n_draws, nrow(fit$m_samples))
 expect_equal(summary(fS)$n_draws, nrow(fS$m_samples))
 expect_false(grepl("noisy", paste(utils::capture.output(print(s)), collapse = "\n")))
 out <- utils::capture.output(postlink:::.print_high_rhat(c(theta = 1.2), n_draws = 100))
 expect_match(out, paste0("Note: split R-hat above 1.05 for theta (1.200); the chain may not have converged within ",
                          "the run, although the split R-hat of fewer than 200 stored draws (here 100) is noisy"),
              fixed = TRUE)
 expect_match(out, "Run more iterations (storing at least 200 draws) and compare fits with several seeds", fixed = TRUE)
 s_thin <- s
 s_thin$n_draws <- 150L
 expect_match(paste(utils::capture.output(print(s_thin)), collapse = "\n"), "(here 150) is noisy", fixed = TRUE)
 # summaries saved without n_draws keep the plain note
 s_old <- s
 s_old$n_draws <- NULL
 expect_match(paste(utils::capture.output(print(s_old)), collapse = "\n"),
              "for theta (1.300); the chain has not converged within the run", fixed = TRUE)

 # fits saved without diagnostics$rhat: the note is computed from the draws,
 # like the Rhat column (here theta drifts during the run)
 old <- fit
 old$diagnostics$rhat <- NULL
 old$diagnostics$ess <- NULL
 old$estimates$theta <- sort(old$estimates$theta)
 so <- summary(old)
 expect_gt(so$theta[, "Rhat"], 1.05)
 expect_equal(so$high_rhat[["theta"]], unname(so$theta[, "Rhat"]))
 expect_true("theta" %in% names(so$low_ess))
 expect_match(paste(utils::capture.output(print(so)), collapse = "\n"), "split R-hat above 1.05 for theta")
})

test_that("a log posterior that is not finite stops the fit or is warned about, also without pilot chains", {
 expect_error(postlink:::.check_lp(c(-Inf, NaN, NA)), "not finite in any stored draw")
 expect_warning(postlink:::.check_lp(c(1, -Inf, 2, 3)), "not finite in 1 of 4 stored draws (25.0%)", fixed = TRUE)
 expect_warning(postlink:::.check_lp(c(1, 2), pilot_lp = c(-Inf, NaN)),
                "any pilot chain.*an outcome or covariates on extreme scales")
 expect_silent(postlink:::.check_lp(c(1, 2), pilot_lp = c(-Inf, 3)))
 expect_silent(postlink:::.check_lp(c(1, 2)))

 # an outcome in extreme units, with pilot chains, without them and with
 # starting values (which skip the pilot chains)
 set.seed(1)
 x <- stats::rnorm(60)
 X <- cbind("(Intercept)" = 1, x = x)
 y <- 1e200 * (1 + x + stats::rnorm(60))
 ctl0 <- list(iterations = 60, burnin.iterations = 30, seed = 1)
 init <- list(beta1 = c(0, 0), beta2 = c(0, 0), theta = 0.5, disp1 = 1, disp2 = 1)
 for (cc in list(ctl0, c(ctl0, pilots = 0), c(ctl0, list(init = init)))) {
  expect_error(suppressWarnings(glmMixBayes(X, y, "gaussian", control = cc)), "not finite in any stored draw")
 }
 # the prior-scale check computed before sampling is reported before the fit stops
 ws <- character()
 expect_error(withCallingHandlers(glmMixBayes(X, y, "gaussian", control = ctl0), warning = function(w) {
  ws <<- c(ws, conditionMessage(w))
  invokeRestart("muffleWarning")
 }), "not finite in any stored draw")
 pw <- grep("prior of", ws, value = TRUE)
 expect_length(pw, 1L)
 expect_match(pw, paste0("^The default prior of '(\\(Intercept\\)|x)' in the correct-match component, ",
                         "normal\\(0, (10|5)\\), does not suit the scale of the data"))
 expect_false(grepl("posterior", pw))
})

# ------------------------------------------------------------------------------
# Acceptance rates and joint moves
# ------------------------------------------------------------------------------
test_that("the acceptance rates of the blocks refer to the sweeps after the burn-in", {
 set.seed(2)
 n <- 200
 x <- stats::rnorm(n)
 yb <- stats::rbinom(n, 1, stats::plogis(0.5 + x))
 X <- cbind("(Intercept)" = 1, x = x)
 Z <- cbind(1, stats::rnorm(n))
 fit_b <- function(iterations, burnin) {
  suppressWarnings(suppressMessages(glmMixBayes(X, yb, "binomial", Z = Z,
                                                control = list(iterations = iterations,
                                                               burnin.iterations = burnin, seed = 1))))
 }
 # one sweep after the burn-in: every rate is 0 or 1
 acc1 <- unlist(fit_b(301, 300)$diagnostics$accept)
 expect_named(acc1, c("beta1", "beta2", "gamma"))
 expect_true(all(acc1 %in% c(0, 1)))
 # ten sweeps after the burn-in: multiples of 1/10
 acc10 <- unlist(fit_b(310, 300)$diagnostics$accept)
 expect_equal(acc10 * 10, round(acc10 * 10))
 # without burn-in every sweep counts: multiples of 1/37
 acc0 <- unlist(fit_b(37, 0)$diagnostics$accept)
 expect_equal(acc0 * 37, round(acc0 * 37))

 # the low-acceptance warning needs at least 50 sweeps in the rates
 expect_silent(postlink:::.check_acceptance(list(beta1 = 0, beta2 = 1, gamma = NA), n_sweeps = 49))
 expect_warning(postlink:::.check_acceptance(list(beta1 = 0, beta2 = 1, gamma = NA), n_sweeps = 50),
                "acceptance rate for the coefficients of component 1")
 # a run with a single sweep after the burn-in has rates of 0 or 1, but no warning
 d <- make_warn(n = 100)
 for (it in c(2, 3)) {
  r <- collect_warnings(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull", adjustment = adjMixBayes(d),
                                  control = list(iterations = it, seed = 1)))
  expect_false(any(grepl("acceptance rate", r$warnings)))
 }
 r <- collect_warnings(glmMixBayes(X, yb, "binomial", Z = Z, control = list(iterations = 301, burnin.iterations = 300,
                                                                          seed = 1)))
 expect_false(any(grepl("acceptance rate", r$warnings)))
})

test_that("summary() notes when the joint moves were off because the burn-in was too short", {
 d <- make_warn()
 X <- cbind("(Intercept)" = 1, x = d$x)
 fit_g <- function(burnin) {
  suppressMessages(glmMixBayes(X, d$y, "gaussian",
                               control = list(iterations = burnin + 150, burnin.iterations = burnin, seed = 1)))
 }
 f40 <- fit_g(40)
 expect_true(all(is.na(unlist(f40$diagnostics$accept_joint))))
 # draws are collected after the first quarter of the burn-in; the proposal
 # needs max(50, 10 P) of them, P = 7 parameters (2 x 2 coefficients, theta,
 # 2 dispersions)
 expect_equal(f40$diagnostics$joint_burnin, c(collected = 30, needed = 70))
 expect_equal(postlink:::.joint_min_burnin(70), 93)
 s <- summary(f40)
 expect_equal(s$joint_short, c(burnin = 40, needed = 93))
 expect_match(paste(utils::capture.output(print(s)), collapse = "\n"),
              "the burn-in of 40 iterations was too short to estimate their proposal; it needs at least 93",
              fixed = TRUE)
 # the smallest sufficient burn-in switches them on
 f92 <- fit_g(92)
 f93 <- fit_g(93)
 expect_true(all(is.na(unlist(f92$diagnostics$accept_joint))))
 expect_false(anyNA(unlist(f93$diagnostics$accept_joint)))
 expect_null(summary(f93)$joint_short)
 # fits saved before joint_burnin was stored: no note
 f40$diagnostics$joint_burnin <- NULL
 expect_null(summary(f40)$joint_short)
 # survival fits
 fS <- suppressMessages(survregMixBayes(X, cbind(d$time, 1), "weibull",
                                        control = list(iterations = 200, burnin.iterations = 20, seed = 1)))
 expect_match(paste(utils::capture.output(print(summary(fS))), collapse = "\n"), "joint moves .* were off")
})

# ------------------------------------------------------------------------------
# Monte Carlo columns of summary()
# ------------------------------------------------------------------------------
test_that("summary() tables add the Monte Carlo standard error, the ESS and the split R-hat", {
 d <- make_warn()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 s <- summary(fit)
 cols <- c("Estimate", "Std. Error", "2.5 %", "97.5 %", "MCSE", "ESS", "Rhat")
 for (nm in c("coefficients", "dispersion", "m.coefficients", "theta")) expect_equal(colnames(s[[nm]]), cols)
 for (nm in c("coefficients2", "dispersion2")) expect_equal(colnames(s[[nm]]), cols[1:4])
 keys <- c("coefficients[(Intercept)]", "coefficients[x]")
 ess <- fit$diagnostics$ess
 expect_equal(unname(s$coefficients[, "ESS"]), unname(round(ess[keys])))
 expect_equal(unname(s$coefficients[, "MCSE"]),
              unname(apply(fit$estimates$coefficients, 2, stats::sd) / sqrt(ess[keys])))
 expect_equal(unname(s$coefficients[, "Rhat"]), unname(fit$diagnostics$rhat[keys]))
 expect_equal(unname(s$theta[, "ESS"]), unname(round(ess["theta"])))
 expect_equal(unname(s$dispersion[, "Rhat"]), unname(fit$diagnostics$rhat["dispersion"]))
 # not stored: computed from the draws
 expect_equal(unname(s$m.coefficients[, "Rhat"]), postlink:::.rhat_split(fit$estimates$m.coefficients[, 1L]))
 out <- paste(utils::capture.output(print(s)), collapse = "\n")
 expect_match(out, "MCSE: Monte Carlo standard error of the estimate; ESS: effective sample size", fixed = TRUE)
 # printed tables: the ESS as a whole number, the split R-hat with three
 # decimals (with significant digits 1.0002 printed as 1), the other columns
 # with significant digits
 tab <- cbind(Estimate = c(a = 1.23456, b = -0.001234), `Std. Error` = c(0.1, 0.2), MCSE = c(0.01, 0.02),
              ESS = c(1234.4, NA), Rhat = c(1.000091, NA))
 pt <- utils::capture.output(postlink:::.print_posterior_table(tab, 4L))
 expect_match(pt[2L], "^a +1\\.234560 +0\\.1 +0\\.01 +1234 +1\\.000$")
 expect_match(pt[3L], "^b +-0\\.001234 +0\\.2 +0\\.02 +NA +NA$")
 tab["a", "Rhat"] <- Inf
 expect_match(utils::capture.output(postlink:::.print_posterior_table(tab, 4L))[2L], "Inf$")
 # tables without the Monte Carlo columns print as before
 tab4 <- tab[, 1:2]
 expect_equal(utils::capture.output(postlink:::.print_posterior_table(tab4, 4L)),
              utils::capture.output(print(tab4, digits = 4L)))
 s1 <- s
 s1$coefficients[, "Rhat"] <- c(1.000091, 1.000226)
 out1 <- utils::capture.output(print(s1))
 expect_true(all(grepl(" 1\\.000$", out1[grep("^(\\(Intercept\\)|x) ", out1)[1:2]])))
 # fits saved without these diagnostics get the same columns from the draws
 old <- fit
 old$diagnostics$ess <- NULL
 old$diagnostics$rhat <- NULL
 so <- summary(old)
 expect_equal(so$coefficients, s$coefficients)
 expect_equal(so$theta, s$theta)

 # linkage covariates and survival models
 fB <- suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl,
                              adjustment = adjMixBayes(d, m.formula = ~ z)))
 sB <- summary(fB)
 expect_equal(colnames(sB$gamma), cols)
 # the split R-hat of gamma equals that of m.coefficients = -gamma
 expect_equal(sB$gamma[, "Rhat"], sB$m.coefficients[, "Rhat"])
 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                  adjustment = adjMixBayes(d), control = ctl))
 sS <- summary(fS, probs = c(0.1, 0.5, 0.9))
 for (nm in c("coef1", "theta", "m.coefficients", "shape1", "inv_shape1")) {
  expect_equal(colnames(sS[[nm]]), c("Estimate", "Std. Error", "10 %", "50 %", "90 %", "MCSE", "ESS", "Rhat"))
 }
 expect_equal(unname(sS$shape1[, "Rhat"]), unname(fS$diagnostics$rhat["shape"]))
 expect_match(paste(utils::capture.output(print(sS)), collapse = "\n"), "Rhat: split R-hat of the single chain",
              fixed = TRUE)
})

test_that("summary() computes the Monte Carlo columns per column when design columns share a name", {
 # cbind(a = 1, a = 2 * x, z) has the column names "a", "a" and "z" (given
 # names are kept, also when duplicated): the two columns named "a" have the
 # same name in diagnostics$ess and diagnostics$rhat, so their values are
 # computed from their own draws
 d <- make_warn()
 fit <- suppressWarnings(suppressMessages(glmMixBayes(cbind(a = 1, a = 2 * d$x, z = d$z), d$y, "gaussian",
                                                      control = ctl)))
 expect_equal(colnames(fit$estimates$coefficients), c("a", "a", "z"))
 b <- fit$estimates$coefficients
 s <- summary(fit)$coefficients
 ess <- apply(b, 2, postlink:::.ess_ar)
 expect_equal(unname(s[, "ESS"]), unname(round(ess)))
 expect_equal(unname(s[, "Rhat"]), unname(apply(b, 2, postlink:::.rhat_split)))
 expect_equal(unname(s[, "MCSE"]), unname(apply(b, 2, stats::sd) / sqrt(ess)))
 # the two unnamed columns differ
 expect_false(isTRUE(all.equal(s[1L, "Rhat"], s[2L, "Rhat"])))
 # stored values are used for uniquely named columns
 expect_equal(s[3L, "Rhat"], fit$diagnostics$rhat[["coefficients[z]"]])
 # and directly: a stored value is used only when its name occurs once
 draws <- cbind(a = stats::rnorm(50), a = stats::rnorm(50), b = stats::rnorm(50))
 mc <- postlink:::.mcmc_columns(draws, "coefficients",
                                ess = c("coefficients[a]" = 1, "coefficients[a]" = 2, "coefficients[b]" = 3),
                                rhat = c("coefficients[a]" = 9, "coefficients[a]" = 9, "coefficients[b]" = 7))
 expect_equal(unname(mc[, "ESS"]), c(round(postlink:::.ess_ar(draws[, 1L])), round(postlink:::.ess_ar(draws[, 2L])), 3))
 expect_equal(unname(mc[, "Rhat"]), c(postlink:::.rhat_split(draws[, 1L]), postlink:::.rhat_split(draws[, 2L]), 7))
})

# ------------------------------------------------------------------------------
# Default burn-in and mi_with()
# ------------------------------------------------------------------------------
test_that("the default burn-in is half of the iterations, at most 1000", {
 b <- function(it) postlink:::.mixbayes_settings(list(), list(iterations = it))$burnin.iterations
 expect_equal(vapply(c(2, 501, 1000, 1001, 1999, 2000, 5000), b, 1L),
              c(1L, 250L, 500L, 500L, 999L, 1000L, 1000L))
 expect_equal(postlink:::.mixbayes_settings(list(iterations = 3000), list())$burnin.iterations, 1000L)
 # an explicit burn-in is kept
 expect_equal(postlink:::.mixbayes_settings(list(), list(iterations = 1500, burnin.iterations = 100))$burnin.iterations,
              100L)
 # the default `control` of the fitting functions holds no burn-in, so the
 # same default applies when `iterations` is given through `...`
 expect_false("burnin.iterations" %in% names(as.list(formals(glmMixBayes)$control)))
 expect_false("burnin.iterations" %in% names(as.list(formals(survregMixBayes)$control)))
 d <- make_warn()
 X <- cbind("(Intercept)" = 1, x = d$x)
 settings <- function(fit) fit$diagnostics$settings[c("iterations", "burnin.iterations")]
 f800 <- suppressWarnings(suppressMessages(glmMixBayes(X, d$y, "gaussian", iterations = 800, seed = 1)))
 expect_equal(settings(f800), list(iterations = 800L, burnin.iterations = 400L))
 expect_equal(nrow(f800$m_samples), 400L)
 f1001 <- suppressWarnings(suppressMessages(glmMixBayes(X, d$y, "gaussian", iterations = 1001, seed = 1)))
 expect_equal(settings(f1001), list(iterations = 1001L, burnin.iterations = 500L))
 fS800 <- suppressWarnings(suppressMessages(survregMixBayes(X, cbind(d$time, d$status), "weibull",
                                                            iterations = 800, seed = 1)))
 expect_equal(settings(fS800), list(iterations = 800L, burnin.iterations = 400L))
 fS1500 <- suppressWarnings(suppressMessages(survregMixBayes(X, cbind(d$time, d$status), "weibull",
                                                             iterations = 1500, seed = 1)))
 expect_equal(settings(fS1500), list(iterations = 1500L, burnin.iterations = 750L))
 # the default of 1e4 iterations keeps a burn-in of 1000
 expect_equal(postlink:::.mixbayes_settings(list(), eval(formals(glmMixBayes)$control))[c("iterations", "burnin.iterations")],
              list(iterations = 10000L, burnin.iterations = 1000L))
})

test_that("mi_with() stops when fewer than two posterior draws can be used", {
 d <- make_warn()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d),
                               control = list(iterations = 2, seed = 1)))
 expect_equal(nrow(fit$m_samples), 1L)
 expect_error(mi_with(fit, data = d), "Only 1 of 1 posterior draws could be used, but at least 2 are needed")
 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                  adjustment = adjMixBayes(d), control = list(iterations = 2, seed = 1)))
 expect_error(mi_with(fS, data = d), "at least 2 are needed")
 expect_error(postlink:::.mi_check_m(1L, 5L, 4L), "or lower `min_n`")
 expect_true(postlink:::.mi_check_m(2L, 5L, 0L))
})
