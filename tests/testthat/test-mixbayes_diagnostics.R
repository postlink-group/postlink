# Diagnostics and conveniences added after the end-to-end review of the
# Bayesian mixture models: effective sample sizes, prior-data conflicts,
# collinear designs, adaptive burn-in, the outcome formula used by mi_with(),
# predict() with named columns, confint() arguments, the identified Weibull
# intercept and the printing of small values.
local_edition(3)

make_diag <- function(n = 150, seed = 3) {
 set.seed(seed)
 x <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, 0.8) == 1
 y <- ifelse(match, 1 + 2 * x + stats::rnorm(n, sd = 0.7), stats::rnorm(n, sd = 2))
 tt <- stats::rweibull(n, shape = 2, scale = exp(ifelse(match, 0.5 + 0.8 * x, 0.2)))
 data.frame(y = y, x = x, x2 = stats::rnorm(n), time = tt, status = 1L)
}
ctl <- list(iterations = 400, burnin.iterations = 200, seed = 4)

test_that("effective sample sizes are stored, warned about for long runs and noted by summary()", {
 d <- make_diag()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 ess <- fit$diagnostics$ess
 expect_named(ess, colnames(posterior_draws(fit)))
 expect_true(all(ess[is.finite(ess)] > 0))
 if (requireNamespace("coda", quietly = TRUE)) {
  # same estimator as coda::effectiveSize()
  expect_equal(unname(ess["coefficients[x]"]),
               unname(coda::effectiveSize(fit$estimates$coefficients[, "x"])), tolerance = 1e-8)
 }

 low <- c("coefficients[(Intercept)]" = 45, "coefficients[x]" = 300, theta = 60, "coefficients2[x]" = 5,
          "m.coefficients[(Intercept)]" = 5)
 expect_warning(postlink:::.check_ess(low, n_draws = 3000), "theta 60")
 expect_silent(postlink:::.check_ess(low, n_draws = 500))
 expect_silent(postlink:::.check_ess(c("coefficients[x]" = 300, "coefficients2[x]" = 5), n_draws = 3000))

 fit$diagnostics$ess <- low
 s <- summary(fit)
 expect_equal(names(s$low_ess), c("coefficients[(Intercept)]", "theta"))
 expect_match(paste(utils::capture.output(print(s)), collapse = "\n"), "low effective sample size")
 fit$diagnostics$ess <- NULL                     # objects fitted with 0.1.2
 expect_length(summary(fit)$low_ess, 0L)

 # with linkage covariates the mismatch model is named m.coefficients[...] in
 # the effective-sample-size check, as in the R-hat check (gamma = -m.coefficients
 # has the same effective sample size and is not listed twice)
 lowB <- c("coefficients[x]" = 300, "m.coefficients[(Intercept)]" = 25, "m.coefficients[z]" = 40,
           "gamma[(Intercept)]" = 25, "gamma[z]" = 40)
 expect_equal(names(postlink:::.low_ess(lowB)), c("m.coefficients[(Intercept)]", "m.coefficients[z]"))
 expect_warning(postlink:::.check_ess(lowB, n_draws = 3000),
                "Low effective sample size: m.coefficients[(Intercept)] 25, m.coefficients[z] 40", fixed = TRUE)
 fB <- suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl,
                              adjustment = adjMixBayes(cbind(d, z = stats::rnorm(nrow(d))), m.formula = ~ z)))
 fB$diagnostics$ess[] <- 500
 fB$diagnostics$ess[c("m.coefficients[z]", "gamma[z]")] <- 30
 sB <- summary(fB)
 expect_equal(names(sB$low_ess), "m.coefficients[z]")
 expect_match(paste(utils::capture.output(print(sB)), collapse = "\n"),
              "low effective sample size for m.coefficients[z] (30)", fixed = TRUE)
 expect_true(all(c("m.coefficients[(Intercept)]", "m.coefficients[z]") %in% names(fB$diagnostics$rhat)))
})

test_that("a coefficient far outside its prior is reported as a prior-data conflict", {
 d <- make_diag()
 d$age <- 60 + 5 * d$y                             # an outcome on the scale of years
 expect_warning(suppressMessages(plglm(age ~ x, family = "gaussian", adjustment = adjMixBayes(d),
                                       control = ctl)),
                "prior standard deviations")
 scaled <- list(intercept1 = "normal(60, 20)", intercept2 = "normal(60, 20)",
                beta1 = "normal(0, 20)", beta2 = "normal(0, 20)")
 expect_no_warning(suppressMessages(plglm(age ~ x, family = "gaussian", control = ctl,
                                          adjustment = adjMixBayes(d, priors = scaled))))
 # an unnamed intercept column (e.g. cbind(1, x)) is called "(Intercept)", other
 # unnamed columns by their position
 expect_warning(suppressMessages(glmMixBayes(cbind(1, d$x), d$age, "gaussian", control = ctl)),
                "default prior of '\\(Intercept\\)'.*default priors are not rescaled to the data")
 # the warning names the help page of the fitting function
 b1 <- cbind("(Intercept)" = rep(50, 10), x = 0)
 ep <- list(beta1_mean = c(0, 0), beta1_sd = c(10, 5))
 expect_warning(postlink:::.check_prior_data_conflict(b1, ep), "in ?glmMixBayes)", fixed = TRUE)
 expect_warning(postlink:::.check_prior_data_conflict(b1, ep, help = "survregMixBayes"),
                "in ?survregMixBayes)", fixed = TRUE)
})

test_that("collinear designs are reported and the default burn-in is half the iterations, at most 1000", {
 d <- make_diag()
 X <- cbind("(Intercept)" = 1, x = d$x, x_copy = 2 * d$x)
 expect_warning(suppressMessages(glmMixBayes(X, d$y, "gaussian", control = ctl)), "x_copy")
 f <- suppressMessages(glmMixBayes(cbind(1, d$x), d$y, "gaussian", control = list(iterations = 500, seed = 1)))
 expect_equal(f$diagnostics$settings$burnin.iterations, 250L)
 expect_equal(nrow(f$m_samples), 250L)
 f2 <- suppressMessages(glmMixBayes(cbind(1, d$x), d$y, "gaussian",
                                    control = list(iterations = 1500, seed = 1)))
 expect_equal(f2$diagnostics$settings$burnin.iterations, 750L)
 expect_error(glmMixBayes(cbind(1, d$x), d$y, "gaussian",
                          control = list(iterations = 500, burnin.iterations = 600)), "smaller")
 expect_error(adjMixBayes(d, priors = list(beta1 = "double_exponential(0, 1)")), "double_exponential")
})

test_that("mi_with() uses the stored outcome formula, not a re-assigned variable", {
 d <- make_diag()
 f <- y ~ x
 fit <- suppressMessages(plglm(f, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))
 f <- y ~ x + x2                                   # the variable is re-used after fitting
 p <- mi_with(fit, data = d)
 expect_equal(names(p$coef), c("(Intercept)", "x"))
 expect_equal(p$coef, mi_with(fit, data = d, formula = y ~ x)$coef)
 rm(f)
 expect_equal(mi_with(fit, data = d)$coef, p$coef)

 fit_fun <- function(dd) {
  g <- y ~ x
  suppressMessages(plglm(g, family = "gaussian", adjustment = adjMixBayes(dd), control = ctl, model = FALSE))
 }
 fm <- fit_fun(d)
 expect_error(mi_with(fm, data = d), "supply `formula`")
 expect_s3_class(mi_with(fm, data = d, formula = y ~ x), "mi_link_pool_glm")
})

test_that("predict() matches named columns of new data and uses log(scale) for Weibull fits", {
 d <- make_diag()
 fit <- suppressMessages(plglm(y ~ x + x2, family = "gaussian", adjustment = adjMixBayes(d), control = ctl,
                               x = TRUE))
 X <- fit$x[1:5, ]
 expect_equal(predict(fit, newx = X[, c("x2", "(Intercept)", "x")]), predict(fit, newx = X))
 expect_error(predict(fit, newx = X[, 1:2]), "2 column")

 fw <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                  adjustment = adjMixBayes(d), control = ctl, x = TRUE))
 b <- fw$estimates$coefficients
 # the reported intercept is the identified (Intercept) + log(scale), so the
 # linear predictor is X beta (log(scale) is not added a second time)
 expect_true(fw$intercept_includes_logscale)
 lp_draws <- fw$x[1:4, ] %*% t(b)
 # named after the rows of the new data, as by predict.glm()
 expect_equal(predict(fw, newdata = fw$x[1:4, ])$component1, rowMeans(lp_draws))
 s <- summary(fw)
 expect_equal(unname(s$coef1["(Intercept)", "Estimate"]), mean(b[, "(Intercept)"]))
 expect_match(paste(utils::capture.output(print(s)), collapse = "\n"),
              "the identified (Intercept) + log(scale)", fixed = TRUE)
 # the identified intercept is close to the truth (0.5), whatever the split
 expect_lt(abs(s$coef1["(Intercept)", "Estimate"] - 0.5), 0.3)
 # fits stored before the intercept included log(scale) keep their predictions
 old <- fw
 old$intercept_includes_logscale <- NULL
 old$estimates$coefficients[, "(Intercept)"] <- b[, "(Intercept)"] - log(fw$estimates$scale)
 old$estimates$coefficients2[, "(Intercept)"] <- fw$estimates$coefficients2[, "(Intercept)"] -
  log(fw$estimates$scale2)
 expect_equal(predict(old, newdata = fw$x[1:4, ]), predict(fw, newdata = fw$x[1:4, ]))
})

test_that("confint.survMixBayes() accepts block = and refuses other arguments", {
 d <- make_diag()
 fw <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ x, dist = "weibull",
                                  adjustment = adjMixBayes(d), control = ctl))
 # block = "<name>" returns that block as a matrix
 expect_equal(confint(fw, block = "theta")["theta", ], confint(fw, parm = "theta")$theta)
 expect_equal(confint(fw, block = "coefficients"), confint(fw)$coef1)
 expect_equal(confint(fw, block = "coef1"), confint(fw)$coef1)
 expect_equal(confint(fw, block = "theta", parm = "theta"), confint(fw, block = "theta"))
 expect_error(confint(fw, block = "nope"), "Unknown block")
 expect_error(confint(fw, type = "x"), "Unused argument")
})

test_that("summaries print small values with significant digits", {
 d <- make_diag()
 d$ysmall <- d$y * 1e-5
 fit <- suppressMessages(plglm(ysmall ~ x, family = "gaussian", control = ctl,
                               adjustment = adjMixBayes(d, priors = list(sigma1 = "cauchy(0, 1e-4)",
                                                                         sigma2 = "cauchy(0, 1e-4)"))))
 out <- paste(utils::capture.output(print(summary(fit))), collapse = "\n")
 expect_match(out, "e-0[45]")
 expect_false(grepl("0\\.0000 +0\\.0000 +0\\.0000", out))
})

test_that("all-safe fits and safe matches that anchor a minority component are reported", {
 d <- make_diag()
 X <- cbind(1, d$x)
 expect_warning(suppressMessages(glmMixBayes(X, d$y, "gaussian", safe.matches = rep(TRUE, nrow(d)), control = ctl)),
                "All records are flagged")
 # three genuine but unusual correct matches flagged as safe (E2E-06): the
 # safe records sit in the tails of the correct-match regression
 set.seed(12)
 n <- 300
 x <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, 0.78) == 1
 y <- ifelse(match, 1 + 2 * x + stats::rnorm(n, sd = 0.5), stats::rnorm(n, sd = 3))
 r <- y - (1 + 2 * x)
 safe <- rep(0L, n)
 safe[which(match)[order(-abs(r[match]))[1:3]]] <- 1L
 ws <- character()
 res <- withCallingHandlers(
  suppressMessages(glmMixBayes(cbind(1, x), y, "gaussian", safe.matches = safe,
                               control = list(iterations = 1500, burnin.iterations = 500, seed = 1))),
  warning = function(w) {
   ws <<- c(ws, conditionMessage(w))
   invokeRestart("muffleWarning")
  })
 # whichever way the chain settles, a component 1 holding a minority of the
 # other records is reported (and only then)
 minority <- mean(res$m_samples[, safe == 0L] == 1L) < 0.5
 expect_identical(any(grepl("atypical", ws)), minority)
 # with a prior saying that correct matches are the minority, no warning
 z_min <- matrix(c(rep(1L, 3), rep(2L, n - 3)), 10, n, byrow = TRUE)
 expect_silent(postlink:::align_mixture_labels(z_min, list(), safe = c(rep(1L, 3), rep(0L, n - 3)),
                                               expect_majority = FALSE))
 expect_warning(postlink:::align_mixture_labels(z_min, list(), safe = c(rep(1L, 3), rep(0L, n - 3)),
                                                expect_majority = TRUE), "atypical")
})

# ------------------------------------------------------------------------------
# Numerical robustness of the C++ sampler (review of the C++ engine)
# ------------------------------------------------------------------------------
test_that("censored gamma survival with a covariate in large units stays finite", {
 set.seed(5)
 N <- 300
 income <- round(stats::rlnorm(N, log(50000), 0.6))
 X <- cbind("(Intercept)" = 1, income = income)
 tt <- stats::rgamma(N, 2, 2 / exp(1 + 1e-5 * income))
 cens <- stats::rexp(N, 1 / 4)
 ev <- as.integer(tt <= cens)
 tt <- pmin(tt, cens)
 f <- suppressWarnings(suppressMessages(survregMixBayes(X, cbind(tt, ev), dist = "gamma",
                                                        control = list(iterations = 1500, burnin.iterations = 500,
                                                                       seed = 2))))
 # the log survival of records whose component-2 mean underflows is -Inf, not NaN
 expect_true(all(is.finite(f$diagnostics$lp)))
 expect_gt(mean(f$estimates$theta), 0.8)                  # no mismatches in these data
 expect_gt(length(unique(f$estimates$theta)), 990)         # theta is updated in every sweep
})

test_that("gaussian fits do not depend on the units of the outcome", {
 set.seed(9)
 N <- 300
 x <- stats::rnorm(N)
 X <- cbind("(Intercept)" = 1, x = x)
 y <- 1 + 2 * x + stats::rnorm(N, sd = 0.5)
 m <- sample(N, 30)
 y[m] <- sample(y[m])
 th <- sapply(c(1, 1e-9), function(k) {
  f <- suppressMessages(glmMixBayes(X, y * k, "gaussian",
                                    control = list(iterations = 1500, burnin.iterations = 500, seed = 1)))
  mean(f$estimates$theta)
 })
 expect_lt(abs(th[1] - th[2]), 0.03)
})

test_that("the C++ ECR relabelling reproduces the R algorithm", {
 ecr_r <- function(z, maxiter = 100L, threshold = 1e-6) {    # the former R implementation
  S <- nrow(z); N <- ncol(z)
  is2 <- z == 2L
  swap <- logical(S)
  cost_prev <- Inf
  for (iter in seq_len(maxiter)) {
   n2 <- colSums(is2[!swap, , drop = FALSE]) + colSums(!is2[swap, , drop = FALSE])
   pivot2 <- n2 > (S - n2)
   n_p2 <- sum(pivot2)
   mism <- rowSums(is2[, !pivot2, drop = FALSE]) + (n_p2 - rowSums(is2[, pivot2, drop = FALSE]))
   new_swap <- mism > (N - mism)
   cost <- sum(pmin(mism, N - mism))
   if (cost > cost_prev) break
   swap <- new_swap
   if (cost_prev - cost <= threshold) break
   cost_prev <- cost
  }
  swap
 }
 set.seed(1)
 for (r in 1:20) {
  S <- sample(2:60, 1); N <- sample(1:40, 1)
  base <- sample(1:2, N, replace = TRUE)
  z <- t(replicate(S, { zz <- base; flip <- stats::runif(N) < 0.2; zz[flip] <- 3L - zz[flip]
                        if (stats::runif(1) < 0.4) 3L - zz else zz }))
  if (N == 1L) z <- matrix(z, ncol = 1L)
  storage.mode(z) <- "integer"
  expect_identical(postlink:::ecr_iterative_two(z), ecr_r(z))
 }
 z1 <- matrix(1L, 5, 4)
 expect_false(any(postlink:::ecr_iterative_two(z1)))
 expect_identical(postlink:::swap_rows_cpp(z1, c(TRUE, FALSE, FALSE, FALSE, TRUE))[c(1, 5), ], matrix(2L, 2, 4))
 expect_equal(postlink:::count_label2_cpp(z1, rep(TRUE, 5), integer(0)), 20)
 expect_equal(postlink:::count_label2_cpp(z1, rep(FALSE, 5), 2:3), 0)
})

test_that("acceptance rates follow the reported component labels", {
 acc <- list(beta1 = 0.1, beta2 = 0.9, gamma = NA)
 expect_equal(postlink:::.relabel_accept(acc, TRUE), list(beta1 = 0.9, beta2 = 0.1, gamma = NA))
 expect_equal(postlink:::.relabel_accept(acc, FALSE), acc)
})

test_that("the joint moves with the indicators integrated out are adapted in the burn-in and reported", {
 d <- make_diag()
 f <- suppressMessages(glmMixBayes(cbind(1, d$x), d$y, "gaussian", control = ctl))
 aj <- f$diagnostics$accept_joint
 expect_named(aj, c("random_walk", "independence"))
 expect_true(all(aj > 0 & aj < 1))
 # without a burn-in the proposals cannot be estimated and the moves stay off
 f0 <- suppressMessages(glmMixBayes(cbind(1, d$x), d$y, "gaussian",
                                    control = list(iterations = 300, burnin.iterations = 0, seed = 1)))
 expect_true(all(is.na(unlist(f0$diagnostics$accept_joint))))
})
