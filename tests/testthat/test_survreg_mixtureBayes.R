local_edition(3)

test_that("survreg_mixtureBayes runs on a simple weibull-like example", {
 set.seed(456)

 # 1) Simulate simple survival data with two mixture components
 n <- 200
 x <- runif(n, -1, 1)
 X <- cbind(1, x)

 beta1 <- c(0,  0.5)
 beta2 <- c(1, -0.5)

 z <- rbinom(n, size = 1, prob = 0.6) + 1

 linpred1 <- as.vector(X %*% beta1)
 linpred2 <- as.vector(X %*% beta2)
 linpred  <- ifelse(z == 1, linpred1, linpred2)

 # Construct positive event times (not necessarily exact Weibull, but valid)
 base_time <- rexp(n, rate = 0.1)
 true_time <- base_time * exp(-linpred)

 # Independent censoring times
 censor_time <- rexp(n, rate = 0.05)
 time   <- pmin(true_time, censor_time)
 status <- as.integer(true_time <= censor_time)  # 1 = event, 0 = censored

 dat <- data.frame(time = time, status = status, x = x)

 # 2) Build survival formula explicitly using survival::Surv
 surv_formula <- survival::Surv(time, status) ~ x

 # 3) Fit the Bayesian survival mixture model (Weibull family)
 #    Suppress Stan warnings (R-hat, ESS) for testing purposes.
 fit <- suppressWarnings(
  survreg_mixtureBayes(
   formula  = surv_formula,
   data     = dat,
   family   = "weibull",
   iterations        = 300,
   burnin.iterations = 150
  )
 )

 # 4) Basic structure checks
 expect_s3_class(fit, "surv_mixtureBayes")
 expect_true(is.matrix(fit$m_samples))
 expect_true(nrow(fit$m_samples) > 0)
 expect_true("estimates" %in% names(fit))
 expect_true(is.matrix(fit$estimates$coefficients))

 # 5) Check that summary runs without error
 expect_error(summary(fit), NA)

 # 6) Check that confint method runs without error
 expect_error(confint.surv_mixtureBayes(fit), NA)

 # 7) Check that the predict method is defined (but do not execute it yet)
 expect_true(is.function(PLDA:::predict.surv_mixtureBayes))
})
