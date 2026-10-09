# Interface and data handling of the Bayesian mixture fits: names of unnamed
# design columns and their use by predict() and summary(), fixed bases of
# data-dependent terms in mi_with(), the rank check of the linkage covariates,
# the checks of the prior specification, the bound on m.rate.sd and the
# printed defaults of adjMixBayes().
local_edition(3)

make_int <- function(n = 160, seed = 11) {
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
# Unnamed design columns
# ------------------------------------------------------------------------------

test_that("unnamed design columns are named once, as '(Intercept)' or 'X<j>'", {
 nm <- function(M, ...) colnames(postlink:::.normalize_design_names(M, ...))
 expect_equal(nm(cbind(1, 2:3)), c("(Intercept)", "X2"))
 expect_equal(nm(unname(cbind(1, 2:3, 4:5))), c("(Intercept)", "X2", "X3"))
 expect_equal(nm(cbind(2:3, 4:5)), c("X1", "X2"))
 # an NA name, a name that is already taken, a named column of ones
 M <- cbind(1, 2:3, 4:5)
 colnames(M) <- c(NA, "X3", "")
 expect_equal(nm(M), c("(Intercept)", "X3", "X3.1"))
 M <- cbind(one = 1, 1, 2:3)
 colnames(M) <- c("(Intercept)", "", "")
 expect_equal(nm(M), c("(Intercept)", "X2", "X3"))
 # names that are given are kept, also when duplicated
 M <- cbind(a = 1:2, a = 3:4)
 expect_equal(nm(M), c("a", "a"))
 expect_equal(nm(cbind(1, 2:3), prefix = "Z"), c("(Intercept)", "Z2"))
})

test_that("glmMixBayes() names unnamed columns and predict() matches by name or position", {
 d <- make_int()
 f <- suppressMessages(glmMixBayes(cbind(1, d$x), d$y, "gaussian", control = ctl))
 expect_equal(colnames(f$estimates$coefficients), c("(Intercept)", "X2"))
 expect_equal(colnames(f$estimates$coefficients2), c("(Intercept)", "X2"))
 expect_true(all(c("coefficients[(Intercept)]", "coefficients[X2]") %in% colnames(posterior_draws(f))))
 expect_true(all(c("coefficients[(Intercept)]", "coefficients[X2]") %in% names(f$diagnostics$ess)))
 expect_equal(rownames(summary(f)$coefficients), c("(Intercept)", "X2"))

 b <- colMeans(f$estimates$coefficients)
 expected <- unname(b[1] + b[2] * c(0, 1))
 # unnamed new data: by position (this failed with 'subscript out of bounds')
 expect_equal(unname(predict(f, newx = cbind(1, c(0, 1)))), expected)
 # reordered columns with the coefficient names: by name
 expect_equal(unname(predict(f, newx = cbind(X2 = c(0, 1), "(Intercept)" = 1))), expected)
 # fits saved with empty coefficient names are matched by position
 fo <- f
 colnames(fo$estimates$coefficients) <- c("", "x")
 expect_equal(unname(predict(fo, newx = cbind(1, c(0, 1)))), expected)
 expect_silent(predict(fo, newx = cbind("(Intercept)" = 1, x = c(0, 1))))
 # by position with a warning when a column bears the name of a coefficient
 # at another position
 expect_warning(p_pos <- predict(fo, newx = cbind(x = c(0, 1), z = 1)),
                "matched to the coefficients by position (the coefficient names are not all non-empty and unique, or not all found in `newx`), but the column \"x\" is at another position",
                fixed = TRUE)
 expect_equal(unname(p_pos), unname(b[1] * c(0, 1) + b[2]))
 expect_warning(predict(fo, newx = cbind(x = c(0, 1), "(Intercept)" = 1)), "check the order of the columns")
 # duplicated coefficient names are matched by position as well
 draws <- matrix(0, 1, 2, dimnames = list(NULL, c("b", "b")))
 newx <- cbind(a = 1:2, b = 3:4)
 expect_equal(postlink:::.align_newx(newx, draws, "newx"), newx)
 expect_error(postlink:::.align_newx(cbind(1:2), draws, "newx"), "has 1 column")

 # unnamed linkage covariates are named "Z<j>"
 fz <- suppressMessages(glmMixBayes(cbind(1, d$x), d$y, "gaussian", control = ctl, Z = cbind(1, d$z)))
 expect_equal(colnames(fz$estimates$gamma), c("(Intercept)", "Z2"))
 expect_equal(names(fz$z_center), "Z2")
})

test_that("survregMixBayes() names unnamed columns and finds the Weibull intercept by its values", {
 d <- make_int()
 y <- cbind(time = d$time, event = d$status)
 fw <- suppressMessages(survregMixBayes(cbind(1, d$x), y, dist = "weibull", control = ctl))
 expect_equal(colnames(fw$estimates$coefficients), c("(Intercept)", "X2"))
 expect_true(fw$intercept_includes_logscale)
 s <- summary(fw)
 expect_equal(rownames(s$coef1), c("(Intercept)", "X2"))
 expect_true(s$intercept_includes_logscale)
 # the identified intercept is reported whatever the name of the column of ones
 X <- cbind(const = 1, x = d$x)
 fc <- suppressMessages(survregMixBayes(X, y, dist = "weibull", control = ctl))
 expect_true(fc$intercept_includes_logscale)
 expect_equal(rownames(summary(fc)$coef1), c("const", "x"))
 # predictions from an unnamed design matrix, by position
 b <- colMeans(fw$estimates$coefficients)
 p <- predict(fw, newdata = cbind(1, c(0, 1)))
 expect_equal(unname(p$component1), unname(b[1] + b[2] * c(0, 1)))
})

# ------------------------------------------------------------------------------
# mi_with(): data-dependent terms keep the basis computed from all records
# ------------------------------------------------------------------------------

make_poly <- function(n = 200, seed = 11) {
 set.seed(seed)
 d <- data.frame(x = stats::rnorm(n, 50, 10))
 d$y <- 2 + 0.3 * (d$x - 50) - 0.02 * (d$x - 50)^2 + stats::rnorm(n)
 mis <- sample(n, 50)
 d$y[mis] <- d$y[sample(mis)]
 lp <- 1 + 0.04 * (d$x - 50) - 0.002 * (d$x - 50)^2
 tt <- stats::rweibull(n, shape = 1.5, scale = exp(lp))
 cc <- stats::rexp(n, 0.05)
 d$time <- pmin(tt, cc)
 d$status <- as.integer(tt <= cc)
 d$time[mis] <- d$time[sample(mis)]
 d
}

test_that("mi_with() refits glm models on the basis of poly() computed from all records", {
 d <- make_poly()
 fit <- collect_warnings(plglm(y ~ poly(x, 2), family = "gaussian", adjustment = adjMixBayes(d),
                               control = ctl))$value
 # the same allocations refitted on the design matrix of all records
 tt <- attr(fit$model, "terms")
 X <- stats::model.matrix(tt, stats::model.frame(tt, d))
 z <- fit$m_samples
 Q <- t(vapply(seq_len(nrow(z)), function(s) {
  i <- z[s, ] == 1L
  stats::lm.fit(X[i, , drop = FALSE], d$y[i])$coefficients
 }, numeric(ncol(X))))

 p_def <- mi_with(fit, data = d)
 expect_equal(names(p_def$coef), colnames(X))
 expect_equal(unname(p_def$coef), unname(colMeans(Q)))
 expect_equal(unname(diag(p_def$B)), unname(apply(Q, 2, stats::var)))
 # an explicit formula gives the basis of all analysed records too
 p_exp <- mi_with(fit, data = d, formula = y ~ poly(x, 2))
 expect_equal(p_exp$coef, p_def$coef)
 # the pooled coefficients are on the scale of the posterior means
 expect_lt(max(abs(p_def$coef - coef(fit)) / abs(coef(fit))), 0.2)
 # glm() refits (a family other than gaussian with identity link) as well
 p_glm <- mi_with(fit, data = d, family = stats::quasi(link = "identity", variance = "constant"))
 expect_equal(p_glm$refit, "glm")
 expect_equal(p_glm$coef, p_def$coef)
 # scale() in a user formula: centred and scaled with all records
 p_sc <- mi_with(fit, data = d, formula = y ~ scale(x))
 Xs <- cbind(1, as.numeric(scale(d$x)))
 Qs <- t(vapply(seq_len(nrow(z)), function(s) {
  i <- z[s, ] == 1L
  stats::lm.fit(Xs[i, , drop = FALSE], d$y[i])$coefficients
 }, numeric(2)))
 expect_equal(unname(p_sc$coef), unname(colMeans(Qs)))
 expect_equal(names(p_sc$coef), c("(Intercept)", "scale(x)"))
})

test_that("mi_with() evaluates data-dependent calls nested in other calls from each draw (as documented)", {
 # R records no prediction call for scale() inside I(): as for predict(), the
 # term is evaluated from the records of each draw, so ?mi_with says to
 # compute such variables in `data` beforehand
 d <- make_poly()
 fit <- collect_warnings(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d), control = ctl))$value
 fit$m_samples[] <- rep(rep(1:2, c(100, nrow(d) - 100)), each = nrow(fit$m_samples))
 p_nested <- suppressWarnings(mi_with(fit, data = d, formula = y ~ I(scale(x)^2)))
 sub <- stats::coef(stats::lm(y ~ I(scale(x)^2), data = d[1:100, ]))
 expect_equal(unname(p_nested$coef), unname(sub))
 # the variable computed beforehand from all records keeps a single basis
 d$sx2 <- as.numeric(scale(d$x))^2
 p_pre <- suppressWarnings(mi_with(fit, data = d, formula = y ~ sx2))
 expect_equal(unname(p_pre$coef), unname(stats::coef(stats::lm(y ~ sx2, data = d[1:100, ]))))
 expect_false(isTRUE(all.equal(unname(p_pre$coef), unname(p_nested$coef))))
})

test_that("mi_with() refits Cox models on the basis of poly() computed from all records", {
 d <- make_poly()
 fit <- collect_warnings(plsurvreg(survival::Surv(time, status) ~ poly(x, 2), dist = "weibull",
                                   adjustment = adjMixBayes(d), control = ctl))$value
 tt <- attr(fit$model, "terms")
 Xf <- stats::model.matrix(tt, stats::model.frame(tt, d))[, -1L]
 z <- fit$m_samples
 Q <- t(vapply(seq_len(nrow(z)), function(s) {
  i <- z[s, ] == 1L
  stats::coef(survival::coxph(survival::Surv(d$time[i], d$status[i]) ~ Xf[i, ], ties = "efron"))
 }, numeric(ncol(Xf))))

 p_def <- mi_with(fit, data = d)
 expect_equal(names(p_def$coef), c("poly(x, 2)1", "poly(x, 2)2"))
 expect_equal(rownames(p_def$vcov), names(p_def$coef))
 expect_equal(unname(p_def$coef), unname(colMeans(Q)), tolerance = 1e-6)
 p_exp <- mi_with(fit, data = d, formula = survival::Surv(time, status) ~ poly(x, 2))
 expect_equal(p_exp$coef, p_def$coef)
 # terms without a data-dependent basis are refitted as written
 p_raw <- mi_with(fit, data = d, formula = survival::Surv(time, status) ~ x + I(x^2))
 expect_equal(names(p_raw$coef), c("x", "I(x^2)"))
})

test_that("the fixed-basis formula of coxph() refits keeps the other variables as written", {
 d <- make_poly()
 d$g <- factor(rep(c("a", "b"), length.out = nrow(d)))
 mf <- stats::model.frame(survival::Surv(time, status) ~ poly(x, 2) * g + I(scale(x)^2), d)
 fx <- postlink:::.mi_fixed_basis(attr(mf, "terms"), d)
 expect_equal(unname(fx$labels), "poly(x, 2)")
 expect_true(".mi_fixed_2_" %in% names(fx$data))
 expect_equal(deparse(fx$formula[[3L]]), ".mi_fixed_2_ * g + I(scale(x)^2)")
 expect_s3_class(fx$formula, "formula")
 nm <- postlink:::.mi_restore_names(c(".mi_fixed_2_1", ".mi_fixed_2_2:gb", "gb"), fx$labels)
 expect_equal(nm, c("poly(x, 2)1", "poly(x, 2)2:gb", "gb"))
 # nothing to fix
 mf <- stats::model.frame(survival::Surv(time, status) ~ x, d)
 expect_length(postlink:::.mi_fixed_basis(attr(mf, "terms"), d)$labels, 0L)
})

test_that("mi_with() keeps penalised pspline() terms in coxph() refits, with the knots of all records", {
 d <- make_poly()
 d$trt <- rep(0:1, length.out = nrow(d))
 # the penalised term stays a call in the formula (coxph() needs it to fit the
 # penalty), with the boundary knots of all records
 mf <- stats::model.frame(survival::Surv(time, status) ~ trt + pspline(x, df = 3), d)
 fx <- postlink:::.mi_fixed_basis(attr(mf, "terms"), d)
 expect_length(fx$labels, 0L)
 ps <- fx$formula[[3L]][[3L]]
 expect_equal(as.character(ps[[1L]]), "pspline")
 expect_equal(eval(ps$Boundary.knots), range(d$x))

 fit <- collect_warnings(plsurvreg(survival::Surv(time, status) ~ trt + x, dist = "weibull",
                                   adjustment = adjMixBayes(d), control = ctl))$value
 z <- fit$m_samples
 # mean of the coxph() refits of formula `f` on the records allocated to the
 # correct matches in each draw
 manual <- function(f, z) {
  K <- length(stats::coef(survival::coxph(f, data = d, ties = "efron")))
  colMeans(t(vapply(seq_len(nrow(z)), function(s) {
   stats::coef(survival::coxph(f, data = d[z[s, ] == 1L, ], ties = "efron"))
  }, numeric(K))))
 }
 bk <- range(d$x)
 p <- mi_with(fit, data = d, formula = survival::Surv(time, status) ~ trt + pspline(x, df = 3))
 expect_equal(names(p$coef), c("trt", paste0("ps(x)", 3:12)))
 expect_equal(unname(p$coef),
              unname(manual(survival::Surv(time, status) ~ trt + pspline(x, df = 3, Boundary.knots = bk), z)),
              tolerance = 1e-6)
 # the penalty arguments are used as written: a fixed theta and method = "aic"
 # (the prediction call of pspline() keeps only the basis arguments and df,
 # so both were refitted with df = 4)
 p_th <- mi_with(fit, data = d, formula = survival::Surv(time, status) ~ trt + pspline(x, theta = 0.9))
 m_th <- manual(survival::Surv(time, status) ~ trt + pspline(x, theta = 0.9, Boundary.knots = bk), z)
 expect_equal(names(p_th$coef), names(m_th))
 expect_equal(p_th$coef, m_th, tolerance = 1e-6)
 p_aic <- mi_with(fit, data = d, formula = survival::Surv(time, status) ~ trt + pspline(x, method = "aic"))
 expect_equal(p_aic$coef,
              manual(survival::Surv(time, status) ~ trt + pspline(x, method = "aic", Boundary.knots = bk), z),
              tolerance = 1e-6)
 # arguments written by position, and survival::pspline(), whose prediction
 # call is not built by makepredictcall()
 mf <- stats::model.frame(survival::Surv(time, status) ~ trt + survival::pspline(x, 3, 0.5), d)
 ps <- postlink:::.mi_fixed_basis(attr(mf, "terms"), d)$formula[[3L]][[3L]]
 expect_equal(ps$df, 3)
 expect_equal(ps$theta, 0.5)
 expect_equal(ps$Boundary.knots, bk)
 # an unpenalised basis is stored as a data column with the knots of all records
 fx <- postlink:::.mi_fixed_basis(
  attr(stats::model.frame(survival::Surv(time, status) ~ pspline(x, df = 3, penalty = FALSE), d), "terms"), d)
 expect_equal(unname(fx$labels), "pspline(x, df = 3, penalty = FALSE)")
 expect_equal(attr(fx$data[[names(fx$labels)]], "Boundary.knots"), bk)

 # the default formula of a fit with a pspline() term (this failed in every
 # draw), with its penalty arguments
 fp <- collect_warnings(plsurvreg(survival::Surv(time, status) ~ pspline(x, theta = 0.9), dist = "weibull",
                                  adjustment = adjMixBayes(d), control = ctl))$value
 p_def <- mi_with(fp, data = d)
 m_def <- manual(survival::Surv(time, status) ~ pspline(x, theta = 0.9, Boundary.knots = bk), fp$m_samples)
 # (nterm = 2.5 * 4 knots with the default df = 4, so 12 coefficients)
 expect_equal(names(p_def$coef), paste0("ps(x)", 3:14))
 expect_equal(p_def$coef, m_def, tolerance = 1e-6)
})

# ------------------------------------------------------------------------------
# Rank check of the linkage covariates
# ------------------------------------------------------------------------------

test_that("linearly dependent linkage covariates are reported", {
 d <- make_int()
 d$z2 <- 2 * d$z
 w <- collect_warnings(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d, m.formula = ~ z + z2),
                             control = ctl))$warnings
 expect_true(any(grepl("match-probability model (Z, from 'm.formula') has linearly dependent columns (z2)",
                       w, fixed = TRUE)))
 # over the records not flagged as safe matches, which inform the model
 d$w <- ifelse(d$safe, stats::rnorm(nrow(d)), 1)
 w <- collect_warnings(plglm(y ~ x, family = "gaussian",
                             adjustment = adjMixBayes(d, m.formula = ~ z + w, safe.matches = safe),
                             control = ctl))$warnings
 expect_true(any(grepl("over the records not flagged as safe matches, which inform it (w)", w, fixed = TRUE)))
 # unnamed columns of Z are named as in the fit
 w <- collect_warnings(glmMixBayes(cbind(1, d$x), d$y, "gaussian", control = ctl,
                                   Z = cbind(1, d$z, 2 * d$z)))$warnings
 expect_true(any(grepl("linearly dependent columns (Z3)", w, fixed = TRUE)))
 # full rank: no warning
 expect_silent(postlink:::.check_Z(cbind(1, d$z), nrow(d)))
 expect_equal(postlink:::.aliased_columns(cbind(1, 1:3, 2 * (1:3))), "column 3")
})

# ------------------------------------------------------------------------------
# Prior specification
# ------------------------------------------------------------------------------

test_that("component-specific prior entries need their component number", {
 expect_null(postlink:::.key_menu("intercept"))
 expect_null(postlink:::.key_menu("beta"))
 expect_null(postlink:::.key_menu("sigma3"))
 expect_null(postlink:::.key_menu("theta1"))
 expect_equal(postlink:::.key_menu("beta2"), "normal")
 expect_equal(postlink:::.key_menu("theta"), "beta")
 expect_equal(postlink:::.key_menu("gamma"), "normal")
 expect_warning(adjMixBayes(NULL, priors = list(intercept = "normal(60, 20)")),
                "need the component number: write intercept1 and/or intercept2.", fixed = TRUE)
 expect_warning(adjMixBayes(NULL, priors = list(beta = "normal(0, 1)", sigma = "cauchy(0, 1)")),
                "write beta1 and/or beta2 (likewise for sigma)", fixed = TRUE)
 expect_warning(postlink:::prepare_mixbayes_priors(list(sigma = "cauchy(0, 1)"), "gaussian", "glm"),
                "write sigma1 and/or sigma2", fixed = TRUE)
 out <- capture_output(suppressWarnings(print(adjMixBayes(NULL, priors = list(beta = "normal(0, 1)",
                                                                              beta1 = "normal(0, 2)")))))
 expect_match(out, "beta\\s+: normal\\(0, 1\\)\\s+\\[not recognised: ignored\\]")
 expect_match(out, "beta1\\s+: normal\\(0, 2\\)\n")
 # the fit does not repeat the warning given when the adjustment object was created
 d <- make_int()
 adj <- suppressWarnings(adjMixBayes(d, priors = list(intercept = "normal(60, 20)", beta1 = "normal(0, 3)")))
 r <- collect_warnings(plglm(y ~ x, family = "gaussian", adjustment = adj, control = ctl))
 expect_false(any(grepl("unrecognised", r$warnings, fixed = TRUE)))
 expect_equal(r$value$priors[["beta1"]], "normal(0, 3)")
 expect_equal(attr(adj$priors, "unrecognised"), "intercept")
 # entries added to the object after its creation are reported at fit time
 adj$priors$beta <- "normal(0, 1)"
 r <- collect_warnings(plglm(y ~ x, family = "gaussian", adjustment = adj, control = ctl))
 expect_true(any(grepl("Ignoring unrecognised prior entries: beta.", r$warnings, fixed = TRUE)))
 expect_false(any(grepl("entries: intercept", r$warnings, fixed = TRUE)))
 adj2 <- adjMixBayes(d)
 adj2$priors <- list(intercept = "normal(60, 20)")
 r <- collect_warnings(plglm(y ~ x, family = "gaussian", adjustment = adj2, control = ctl))
 expect_true(any(grepl("Ignoring unrecognised prior entries: intercept.", r$warnings, fixed = TRUE)))
 expect_equal(r$value$priors[["intercept1"]], "normal(0, 10)")
 # priors given to plglm() itself are checked when the model is fitted
 r <- collect_warnings(plglm(y ~ x, family = "gaussian", adjustment = adjMixBayes(d),
                             priors = list(intercept = "normal(0, 1)"), control = ctl))
 expect_true(any(grepl("Ignoring unrecognised prior entries: intercept.", r$warnings, fixed = TRUE)))
})

test_that("intercept priors without an intercept column are reported and left out of fit$priors", {
 d <- make_int()
 adj <- adjMixBayes(d, priors = list(intercept1 = "normal(1, 1)", intercept2 = "normal(0, 1)",
                                     beta1 = "normal(0, 3)"))
 r <- collect_warnings(plglm(y ~ 0 + x, family = "gaussian", adjustment = adj, control = ctl))
 expect_true(any(grepl(paste0("The intercept1 and intercept2 priors are not used: the design matrix has no ",
                              "intercept column"), r$warnings, fixed = TRUE)))
 expect_false(any(c("intercept1", "intercept2") %in% names(r$value$priors)))
 expect_equal(r$value$priors[["beta1"]], "normal(0, 3)")
 # with an intercept column they are used and reported
 fI <- collect_warnings(plglm(y ~ x, family = "gaussian", adjustment = adj, control = ctl))
 expect_false(any(grepl("not used", fI$warnings)))
 expect_equal(fI$value$priors[["intercept1"]], "normal(1, 1)")
 # survival fits
 rS <- collect_warnings(survregMixBayes(cbind(x = d$x), cbind(time = d$time, event = d$status),
                                        dist = "weibull", priors = list(intercept1 = "normal(0, 1)"),
                                        control = ctl))
 expect_true(any(grepl("The intercept1 prior is not used", rS$warnings, fixed = TRUE)))
 expect_false("intercept1" %in% names(rS$value$priors))
 # a column of ones that is not first is named "(Intercept)" but receives the
 # beta1 / beta2 prior, which the warning says
 r1 <- collect_warnings(glmMixBayes(unname(cbind(d$x, 1)), d$y, "gaussian",
                                    priors = list(intercept1 = "normal(1, 1)"), control = ctl))
 expect_equal(colnames(r1$value$estimates$coefficients), c("X1", "(Intercept)"))
 expect_true(any(grepl(paste0("The intercept1 prior is not used: the intercept prior applies only to a first ",
                              "column of ones, and the column of ones of the design matrix (column 2, ",
                              "\"(Intercept)\") is not first"), r1$warnings, fixed = TRUE)))
 expect_false(any(grepl("has no intercept column", r1$warnings, fixed = TRUE)))
 expect_false("intercept1" %in% names(r1$value$priors))
 # helper: nothing to drop
 expect_identical(postlink:::.drop_unused_intercept_priors(list(intercept1 = "normal(0, 1)"), TRUE),
                  list(intercept1 = "normal(0, 1)"))
 expect_null(postlink:::.drop_unused_intercept_priors(NULL, FALSE))
})

test_that("the precedence rules among the match-probability priors are reported", {
 expect_message(postlink:::prepare_mixbayes_priors(list(theta = "beta(2, 2)", gamma_intercept = "normal(1, 1)"),
                                                   "gaussian", "glm", use_logistic = TRUE),
                "The 'theta' prior is not used")
 expect_message(pr <- postlink:::prepare_mixbayes_priors(list(gamma = "normal(0, 1)", gamma_slope = "normal(0, 2)"),
                                                         "gaussian", "glm", use_logistic = TRUE),
                "'gamma' prior (an alias of 'gamma_slope') is not used", fixed = TRUE)
 expect_equal(pr$prior_gamma_slope_sd, 2)
 # no message when nothing is overridden
 expect_silent(postlink:::prepare_mixbayes_priors(list(theta = "beta(2, 2)"), "gaussian", "glm",
                                                  use_logistic = TRUE))
 expect_silent(postlink:::prepare_mixbayes_priors(list(gamma = "normal(0, 1)"), "gaussian", "glm",
                                                  use_logistic = TRUE))
 out <- capture_output(print(adjMixBayes(NULL, m.formula = ~ z,
                                         priors = list(theta = "beta(2, 2)", gamma_intercept = "normal(1, 1)",
                                                       gamma = "normal(0, 1)", gamma_slope = "normal(0, 2)"))))
 expect_match(out, "user-specified (the theta prior is not used)", fixed = TRUE)
 expect_match(out, "its alias gamma is not used", fixed = TRUE)
})

test_that("the (m.rate, m.rate.sd) pair is checked only when m.rate is used", {
 expect_error(adjMixBayes(NULL, m.rate = 0.01), "strictly less than")
 a <- adjMixBayes(NULL, m.rate = 0.01, priors = list(theta = "beta(2, 2)"))
 expect_equal(a$m.rate, 0.01)
 expect_silent(adjMixBayes(NULL, m.formula = ~ z, m.rate = 0.01,
                           priors = list(gamma_intercept = "normal(2, 1)")))
 # without linkage covariates gamma_intercept does not replace m.rate
 expect_error(suppressWarnings(adjMixBayes(NULL, m.rate = 0.01, priors = list(gamma_intercept = "normal(2, 1)"))),
              "strictly less than")
 # each value is still checked on its own
 expect_error(adjMixBayes(NULL, m.rate = 1.5, priors = list(theta = "beta(2, 2)")), "strictly between")
 expect_error(adjMixBayes(NULL, m.rate = 0.01, m.rate.sd = -1, priors = list(theta = "beta(2, 2)")), "positive")
 expect_message(pr <- postlink:::prepare_mixbayes_priors(list(theta = "beta(2, 2)"), "gaussian", "glm",
                                                         m.rate = 0.01), "not used")
 expect_equal(c(pr$prior_theta_alpha, pr$prior_theta_beta), c(2, 2))
 expect_error(postlink:::prepare_mixbayes_priors(NULL, "gaussian", "glm", m.rate = 0.01), "strictly less than")
 d <- make_int()
 fit <- suppressMessages(plglm(y ~ x, family = "gaussian", control = ctl,
                               adjustment = adjMixBayes(d, m.rate = 0.01, priors = list(theta = "beta(2, 2)"))))
 expect_equal(fit$priors[["theta"]], "beta(2, 2)")
})

test_that("fit$priors keeps the spelling of exponential priors", {
 out <- postlink:::prepare_mixbayes_priors(list(phi1 = "exponential(3)"), "gamma", "glm")
 expect_equal(out$prior_phi1_dist, "gamma")
 expect_equal(out$prior_phi1_args, c(1, 3))
 s <- postlink:::.prior_strings(out, "gamma", "glm", FALSE)
 expect_equal(s[["phi1"]], "exponential(3)")
 expect_equal(s[["phi2"]], "gamma(2, 0.1)")
 # the survival gamma defaults
 def <- postlink:::prepare_mixbayes_priors(NULL, "gamma", "survival")
 expect_equal(postlink:::.prior_strings(def, "gamma", "survival", FALSE)[["phi1"]], "exponential(1)")
 # gamma(1, 1) is the same distribution as the default exponential(1)
 pf <- postlink:::prepare_mixbayes_priors(list(phi1 = "gamma(1, 1)"), "gamma", "survival")
 expect_true(postlink:::.default_component_priors(pf, "gamma", "survival"))
 expect_equal(postlink:::.prior_strings(pf, "gamma", "survival", FALSE)[["phi1"]], "gamma(1, 1)")
 # without an intercept column, intercept1 / intercept2 are left out
 expect_named(postlink:::.prior_strings(def, "gamma", "survival", FALSE, intercept = FALSE),
              c("beta1", "beta2", "phi1", "phi2", "theta"))
})

# ------------------------------------------------------------------------------
# The bound on m.rate.sd and the printed defaults
# ------------------------------------------------------------------------------

test_that("the m.rate.sd bound is rounded down and the warnings describe the prior used", {
 expect_equal(postlink:::.mrate_sd_bound_text(0.05), "0.04755")
 expect_equal(postlink:::.floor_signif(0.0475594866, 4), 0.04755)
 expect_equal(postlink:::.floor_signif(0.25, 4), 0.25)
 for (m in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.95)) {
  b <- as.numeric(postlink:::.mrate_sd_bound_text(m))
  expect_lte(b, postlink:::.mrate_sd_unimodal(m))
  bp <- postlink:::mrate_to_beta(m, b)
  expect_gt(min(bp$alpha, bp$beta), 1)
 }
 # without linkage covariates: the beta prior on theta
 expect_warning(postlink:::prepare_mixbayes_priors(NULL, "gaussian", "glm", m.rate = 0.05),
                paste0("implies a beta(3.56, 0.187) prior on the match probability, whose density is ",
                       "unbounded at 1 because a shape parameter is below 1"), fixed = TRUE)
 expect_warning(postlink:::prepare_mixbayes_priors(NULL, "gaussian", "glm", m.rate = 0.05),
                "With m.rate.sd at most 0.04755 both shape parameters exceed 1", fixed = TRUE)
 # with linkage covariates: the normal prior on the logit scale actually used
 w <- tryCatch(postlink:::prepare_mixbayes_priors(NULL, "gaussian", "glm", use_logistic = TRUE, m.rate = 0.05),
               warning = conditionMessage)
 expect_match(w, "the prior normal(6.76, 5.48) on the logit scale", fixed = TRUE)
 expect_match(w, "prior median 0.999, 95% interval (0.0183, 1) and standard deviation 0.267", fixed = TRUE)
 expect_match(w, "(rather than m.rate.sd = 0.1)", fixed = TRUE)
 expect_match(w, "With m.rate.sd = 0.04755 this standard deviation would be 0.0767", fixed = TRUE)
 expect_false(grepl("beta(", w, fixed = TRUE))
 # no warning at the printed bound, nor at the exact bound, where the smaller
 # shape parameter is 1 up to rounding error
 expect_silent(postlink:::prepare_mixbayes_priors(NULL, "gaussian", "glm", m.rate = 0.05, m.rate.sd = 0.04755))
 for (m in c(0.05, 0.2, 0.7)) {
  sd_m <- postlink:::.mrate_sd_unimodal(m)
  expect_silent(postlink:::prepare_mixbayes_priors(NULL, "gaussian", "glm", m.rate = m, m.rate.sd = sd_m))
  expect_silent(postlink:::prepare_mixbayes_priors(NULL, "gaussian", "glm", use_logistic = TRUE,
                                                   m.rate = m, m.rate.sd = sd_m))
  expect_false(grepl("unbounded", capture_output(print(adjMixBayes(NULL, m.rate = m, m.rate.sd = sd_m)))))
 }
 expect_equal(postlink:::.shape_below_one(c(1 - 1e-12, 0.99, 1.2)), c(FALSE, TRUE, FALSE))
 # the notes of print.adjMixBayes()
 out <- capture_output(print(adjMixBayes(NULL, m.rate = 0.05)))
 expect_match(out, "density is unbounded")
 expect_match(out, "at 1; with m.rate.sd at most 0.04755 both shape parameters exceed 1.", fixed = TRUE)
 out <- capture_output(print(adjMixBayes(NULL, m.formula = ~ z, m.rate = 0.05)))
 expect_match(out, "on the probability scale this prior has standard deviation 0.267")
 expect_match(out, "with m.rate.sd = 0.04755 it would be 0.0767", fixed = TRUE)
 # the probability-scale moments of a normal prior on the logit scale
 mom <- postlink:::.logitnormal_moments(1.5, 0.3)
 set.seed(1)
 p <- stats::plogis(stats::rnorm(2e5, 1.5, 0.3))
 expect_equal(unname(mom), c(mean(p), stats::sd(p)), tolerance = 0.01)
})

test_that("print.adjMixBayes() says that the listed defaults are used", {
 out <- capture_output(print(adjMixBayes(NULL)))
 expect_match(out, "None specified; the defaults below are used.", fixed = TRUE)
 expect_false(grepl("symmetric", out))
 expect_match(out, "phi1, phi2 ~ exponential(1)", fixed = TRUE)
})
