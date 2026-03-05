local_edition(3)

# Helper function to simulate GLM data with ELE linkage errors
# Adapted from Chambers (2009) simulation logic
generate_ele_glm_data <- function(n = 200, m_rate = 0.1, family = "gaussian") {
 set.seed(42)

 # 1. Design Matrix and True Parameters
 x1 <- rnorm(n)
 x2 <- rbinom(n, 1, 0.5)
 X <- cbind(Intercept = 1, x1 = x1, x2 = x2)
 true_beta <- c(1.5, 0.8, -0.5)

 # 2. Generate Response based on family
 eta <- as.vector(X %*% true_beta)

 if (family == "gaussian") {
  y <- eta + rnorm(n, sd = 0.5)
  fam_obj <- gaussian()
 } else if (family == "binomial") {
  prob <- plogis(eta)
  y <- rbinom(n, 1, prob)
  fam_obj <- binomial()
 } else if (family == "poisson") {
  mu <- exp(eta)
  y <- rpois(n, mu)
  fam_obj <- poisson()
 }

 # 3. Induce Linkage Errors (ELE model within a single block)
 # Shuffling a portion of the response vector to simulate mismatches
 mis_idx <- sample(1:n, size = round(m_rate * n))
 linked_y <- y
 if (length(mis_idx) > 1) {
  linked_y[mis_idx] <- y[sample(mis_idx)]
 }

 return(list(
  X = X,
  y = linked_y,
  true_beta = true_beta,
  family = family,
  fam_obj = fam_obj
 ))
}

test_that("glmELE engine fits all weighting methods for Gaussian", {
 dat <- generate_ele_glm_data(n = 300, m_rate = 0.15, family = "gaussian")

 # Test with "all" weight matrices
 fit <- glmELE(
  x = dat$X,
  y = dat$y,
  family = "gaussian",
  m.rate = 0.15,
  blocks = rep(1, 300),
  weight.matrix = "all"
 )

 expect_s3_class(fit, "glmELE")
 expect_equal(nrow(fit$coefficients), 3) # ratio, LL, BLUE
 expect_true(all(c("ratio", "LL", "BLUE") %in% rownames(fit$coefficients)))
 expect_true(all(diag(fit$var$BLUE) > 0)) # Check variance matrix
})

test_that("glmELE works for Binomial and Poisson families", {
 # Binomial
 dat_b <- generate_ele_glm_data(n = 200, m_rate = 0.1, family = "binomial")
 fit_b <- glmELE(dat_b$X, dat_b$y, family = "binomial", m.rate = 0.1, blocks = rep(1, 200))
 expect_equal(fit_b$family$family, "binomial")

 # Poisson
 dat_p <- generate_ele_glm_data(n = 200, m_rate = 0.1, family = "poisson")
 fit_p <- glmELE(dat_p$X, dat_p$y, family = "poisson", m.rate = 0.1, blocks = rep(1, 200))
 expect_equal(fit_p$family$family, "poisson")
})

test_that("S3 methods (print, summary, vcov, confint) behave as expected", {
 dat <- generate_ele_glm_data(n = 150)

 adj <- adjELE(linked.data = data.frame(dat$y, dat$X),
               m.rate = 0.1,
               weight.matrix = "all")

 fit <- plglm(dat.y ~ x1 + x2, data = data.frame(dat.y = dat$y, dat$X),
              family = "gaussian", adjustment = adj)

 expect_output(print(fit))
 s_fit <- summary(fit)
 expect_s3_class(s_fit, "summary.glmELE")

 expect_true("BLUE" %in% names(s_fit$coefficients))

 expect_equal(dim(vcov(fit, weight.matrix = "BLUE")), c(3, 3))

 ci <- confint(fit, level = 0.95, weight.matrix = "LL")
 expect_equal(nrow(ci), 3)
 expect_equal(colnames(ci), c("2.5 %", "97.5 %"))
})

test_that("Predict method handles new data and standard errors", {
 dat <- generate_ele_glm_data(n = 100)
 adj <- adjELE(linked.data = data.frame(dat.y = dat$y, x1 = dat$X[,2], x2 = dat$X[,3]), m.rate = 0.1)
 fit <- plglm(dat.y ~ x1 + x2, data = data.frame(dat.y = dat$y, x1 = dat$X[,2], x2 = dat$X[,3]),
              family = "gaussian", adjustment = adj)

 # In-sample prediction
 pred_in <- predict(fit, type = "response")
 expect_length(pred_in, 100)

 # New data with SE and intervals
 new_df <- data.frame(x1 = c(0, 0.5), x2 = c(1, 0))
 pred_new <- predict(fit, newdata = new_df, se.fit = TRUE, interval = "confidence")

 expect_equal(nrow(pred_new$fit), 2)
 expect_equal(colnames(pred_new$fit), c("fit", "lwr", "upr"))
 expect_length(pred_new$se.fit, 2)
 expect_true(pred_new$residual.scale > 0)
})

test_that("glmELE handles multi-block m.rates and audit.size", {
 # Simulate 2 blocks with different rates
 n <- 200
 blocks <- c(rep(1, 100), rep(2, 100))
 dat <- generate_ele_glm_data(n = n)

 # m.rate as a vector
 fit <- glmELE(
  x = dat$X, y = dat$y,
  m.rate = c(0.1, 0.05),
  blocks = blocks,
  audit.size = c(20, 20) # Variance inflation from estimated alpha
 )

 expect_equal(nrow(fit$coefficients), 3)
 expect_true(!is.null(fit$var$BLUE))
})

test_that("fitglm.adjELE correctly synchronizes subsets and row names", {
 n <- 100
 dat <- generate_ele_glm_data(n = n)
 df <- data.frame(y = dat$y, x1 = dat$X[,2], x2 = dat$X[,3])
 rownames(df) <- paste0("ID_", 1:n)

 adj <- adjELE(linked.data = df, m.rate = 0.1, weight.matrix = "all")

 # Fit on a subset
 # fitglm.adjELE must find row names ID_11 to ID_60 in the linked.data
 fit_sub <- plglm(y ~ x1 + x2, data = df, subset = 11:60, adjustment = adj)

 expect_equal(fit_sub$df.residual, 50 - 3)
 expect_equal(length(fit_sub$residuals["BLUE",]), 50)
 expect_equal(names(fit_sub$residuals["BLUE",])[1], "ID_11")
})

test_that("Error handling for invalid inputs", {
 dat <- generate_ele_glm_data(n = 20)

 # Dimension mismatch
 expect_error(glmELE(x = dat$X[1:10,], y = dat$y, m.rate = 0.1, blocks = rep(1, 20)))

 # Binomial response not 0/1
 expect_error(glmELE(x = dat$X, y = runif(20, 2, 5), family = "binomial", m.rate = 0.1, blocks = rep(1, 20)))

 # Unsupported family
 expect_error(glmELE(x = dat$X, y = dat$y, family = "inverse.gaussian", m.rate = 0.1, blocks = rep(1, 20)))
})
