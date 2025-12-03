local_edition(3)

test_that("glm_mixtureBayes runs on a simple gaussian example", {
 set.seed(123)

 # 1) Simulate a simple two-component Gaussian mixture with one covariate
 n <- 200
 x <- runif(n, -2, 2)
 X <- cbind(1, x)  # intercept + one covariate

 beta1 <- c(0, 1)   # component 1 coefficients
 beta2 <- c(2, -1)  # component 2 coefficients

 z <- rbinom(n, size = 1, prob = 0.6) + 1  # labels in {1, 2}
 mu1 <- as.vector(X %*% beta1)
 mu2 <- as.vector(X %*% beta2)
 mu  <- ifelse(z == 1, mu1, mu2)

 y <- rnorm(n, mean = mu, sd = 1)
 dat <- data.frame(y = y, x = x)

 # 2) Fit the Bayesian GLM mixture model (Gaussian family)
 #    Suppress Stan warnings about ESS etc. for testing purposes.
 fit <- suppressWarnings(
  glm_mixtureBayes(
   formula = y ~ x,
   data    = dat,
   family  = "gaussian",
   iterations        = 300,   # small numbers for testing speed
   burnin.iterations = 150
  )
 )

 # 3) Basic structure checks
 expect_s3_class(fit, "glm_mixtureBayes")
 expect_true(is.matrix(fit$m_samples))
 expect_true(nrow(fit$m_samples) > 0)
 expect_true("estimates" %in% names(fit))
 expect_true(is.matrix(fit$estimates$coefficients))

 # 4) Check that summary runs without error
 expect_error(summary(fit), NA)

 # 5) Check that confint method runs without error
 expect_error(confint.glm_mixtureBayes(fit), NA)

 # 6) Check that the predict method is defined (but do not execute it yet)
 expect_true(is.function(PLDA:::predict.glm_mixtureBayes))
})
