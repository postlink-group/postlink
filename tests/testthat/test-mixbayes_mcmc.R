# Tests for posterior_draws(), the coda::as.mcmc() methods and intercept-only
# designs of the Bayesian mixture models.
local_edition(3)

make_mcmc_data <- function(n = 120, seed = 3) {
 set.seed(seed)
 x <- stats::rnorm(n)
 z <- stats::rnorm(n)
 match <- stats::rbinom(n, 1, 0.8) == 1
 y <- ifelse(match, 1 + 2 * x + stats::rnorm(n, sd = 0.5), stats::rnorm(n, sd = 2))
 tt <- stats::rweibull(n, 2, exp(ifelse(match, 0.5 + 0.8 * x, 0.2)))
 cens <- stats::rexp(n, 0.1)
 list(
  X = cbind("(Intercept)" = 1, x = x),
  Z = cbind("(Intercept)" = 1, z = z),
  y = y,
  time = pmin(tt, cens),
  status = as.integer(tt <= cens)
 )
}

ctl_thin <- list(iterations = 400, burnin.iterations = 200, thin = 4, seed = 5)

test_that("posterior_draws() collects every stored block in chain order", {
 d <- make_mcmc_data()
 fit <- suppressMessages(glmMixBayes(d$X, d$y, "gaussian", control = ctl_thin))
 pd <- posterior_draws(fit)
 expect_equal(dim(pd), c(50L, 9L))
 expect_equal(colnames(pd), c("coefficients[(Intercept)]", "coefficients[x]",
                              "coefficients2[(Intercept)]", "coefficients2[x]",
                              "m.coefficients[(Intercept)]",
                              "theta", "dispersion", "dispersion2", "lp"))
 expect_equal(unname(pd[, "theta"]), fit$estimates$theta)
 expect_equal(unname(pd[, "lp"]), fit$diagnostics$lp)
 expect_equal(unname(pd[, 1:2]), unname(fit$estimates$coefficients))

 # unnamed design matrix: the columns are named "(Intercept)" (a column of
 # ones) and "X<j>" by glmMixBayes()
 fu <- suppressMessages(glmMixBayes(unname(d$X), d$y, "gaussian", control = ctl_thin))
 expect_equal(colnames(posterior_draws(fu))[1:2], c("coefficients[(Intercept)]", "coefficients[X2]"))

 # logistic match-probability model: gamma block instead of theta
 fB <- suppressMessages(glmMixBayes(d$X, d$y, "gaussian", Z = d$Z, control = ctl_thin))
 pB <- posterior_draws(fB)
 expect_true(all(c("gamma[(Intercept)]", "gamma[z]") %in% colnames(pB)))
 expect_false("theta" %in% colnames(pB))
 expect_equal(unname(pB[, c("m.coefficients[(Intercept)]", "m.coefficients[z]")]),
              unname(-pB[, c("gamma[(Intercept)]", "gamma[z]")]))

 # survival: shape and (Weibull) scale blocks
 fW <- suppressMessages(survregMixBayes(d$X, survival::Surv(d$time, d$status), "weibull",
                                        control = ctl_thin))
 expect_equal(colnames(posterior_draws(fW)),
              c("coefficients[(Intercept)]", "coefficients[x]",
                "coefficients2[(Intercept)]", "coefficients2[x]",
                "m.coefficients[(Intercept)]",
                "theta", "shape", "shape2", "scale", "scale2", "lp"))
})

test_that("coda::as.mcmc() returns the draws with the MCMC settings as attributes", {
 skip_if_not_installed("coda")
 d <- make_mcmc_data()
 fit <- suppressMessages(glmMixBayes(d$X, d$y, "gaussian", control = ctl_thin))
 m <- coda::as.mcmc(fit)
 expect_s3_class(m, "mcmc")
 expect_equal(attr(m, "mcpar"), c(204, 400, 4))
 expect_equal(as.numeric(m[, "theta"]), fit$estimates$theta)
 expect_true(all(is.finite(coda::effectiveSize(m))))

 fS <- suppressMessages(survregMixBayes(d$X, survival::Surv(d$time, d$status), "gamma",
                                        Z = d$Z, control = ctl_thin))
 expect_equal(attr(coda::as.mcmc(fS), "mcpar"), c(204, 400, 4))
})

test_that("the comparison of several seeds recommended in ?glmMixBayes runs with linkage covariates", {
 skip_if_not_installed("coda")
 # gamma and m.coefficients = -gamma are both in the draws, so their
 # covariance matrix is singular: the help page uses multivariate = FALSE
 d <- make_mcmc_data()
 ctl2 <- function(seed) list(iterations = 400, burnin.iterations = 200, seed = seed)
 f1 <- suppressWarnings(suppressMessages(glmMixBayes(d$X, d$y, "gaussian", Z = d$Z, control = ctl2(1))))
 f2 <- suppressWarnings(suppressMessages(glmMixBayes(d$X, d$y, "gaussian", Z = d$Z, control = ctl2(2))))
 g <- coda::gelman.diag(coda::mcmc.list(coda::as.mcmc(f1), coda::as.mcmc(f2)), multivariate = FALSE)
 expect_equal(rownames(g$psrf), colnames(posterior_draws(f1)))
 expect_true(all(is.finite(g$psrf[, "Point est."])))
 expect_null(g$mpsrf)
 # the same for survival fits
 s1 <- suppressWarnings(suppressMessages(survregMixBayes(d$X, survival::Surv(d$time, d$status), "weibull",
                                                         Z = d$Z, control = ctl2(1))))
 s2 <- suppressWarnings(suppressMessages(survregMixBayes(d$X, survival::Surv(d$time, d$status), "weibull",
                                                         Z = d$Z, control = ctl2(2))))
 gs <- coda::gelman.diag(coda::mcmc.list(coda::as.mcmc(s1), coda::as.mcmc(s2)), multivariate = FALSE)
 expect_equal(rownames(gs$psrf), colnames(posterior_draws(s1)))
})

test_that("intercept-only designs keep their matrix structure through every method", {
 d <- make_mcmc_data()
 dat <- data.frame(y = d$y, time = d$time, status = d$status)
 fit <- suppressMessages(plglm(y ~ 1, family = "gaussian", adjustment = adjMixBayes(dat),
                               control = ctl_thin))
 expect_equal(dim(fit$estimates$coefficients), c(50L, 1L))
 expect_equal(colnames(fit$estimates$coefficients), "(Intercept)")
 expect_equal(rownames(summary(fit)$coefficients), "(Intercept)")
 expect_equal(dim(confint(fit)), c(1L, 2L))
 expect_equal(dim(vcov(fit)), c(1L, 1L))
 expect_length(predict(fit, newx = matrix(1, 2, 1)), 2L)
 expect_named(coef(fit), "(Intercept)")

 fS <- suppressMessages(plsurvreg(survival::Surv(time, status) ~ 1, dist = "weibull",
                                  adjustment = adjMixBayes(dat), control = ctl_thin))
 expect_equal(dim(fS$estimates$coefficients), c(50L, 1L))
 expect_equal(dim(confint(fS)$coef1), c(1L, 2L))
 expect_equal(dim(summary(fS)$coef1), c(1L, 7L))
 expect_length(predict(fS, newdata = matrix(1, 2, 1))$component1, 2L)
 expect_output(print(fS), "\\(Intercept\\)")
})
