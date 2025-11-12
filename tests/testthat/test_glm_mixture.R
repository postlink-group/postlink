local_edition(3)

test_that("glm_mixture example 1 works", {
  
  ################################################################################
  ### 1. Generate the Test Data 
  ################################################################################
  set.seed(123)
  n <- 1000
  
  x <- seq(from = -3, to = 3, length = n) # continuous - uniformly b/w [-3,3]
  d <- c(rep(0, n/2), rep(1, n/2)) # binary covariate - 0 or 1 equally
  X <- cbind(1, d,x, x*d) # design matrix of covariates 
  
  Delta <- matrix(nrow = n, ncol = 1, data = 1) # intercept-only model
  betastar <- c(0.5, -1.5, 1, 0.5)
  mustar <- X %*% betastar
  sigma <- 0.25
  alphastar <- 0.05
  
  y <- rnorm(n = n, mean = mustar, sd = sigma)
  m <- rbinom(n, size = 1, prob = alphastar)
  shuffled_ix <- sample(which(m == 1))
  yperm <- y
  yperm[shuffled_ix] <- yperm[c(shuffled_ix[2:length(shuffled_ix)], shuffled_ix[1])]
  
  ################################################################################
  ### 2. Perform adjustment method
  ################################################################################

  expect_snapshot(glm_mixture(formula = yperm ~ X - 1))
  
  fit <- glm_mixture(formula = yperm ~ X - 1) 
  expect_snapshot(attributes(fit))
  expect_snapshot(print(fit))
  expect_snapshot(summary(fit))
  expect_snapshot(vcov(fit))
  expect_snapshot(confint(fit))
  expect_snapshot(predict(fit, se.fit = TRUE, interval = "confidence"))
  
})