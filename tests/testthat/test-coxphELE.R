local_edition(3)

# Helper function to simulate CoxPH data with ELE linkage errors
# Based on Vo et al. (2024) logic
generate_ele_cox_data <- function(n = 200, m_rate = 0.15, num_blocks = 1) {
 set.seed(123)

 # 1. Simulate Covariates and Survival Times
 x1 <- rnorm(n)
 x2 <- rbinom(n, 1, 0.5)
 true_beta <- c(x1 = 0.8, x2 = -0.5)

 true_hazard <- exp(0.8 * x1 - 0.5 * x2)
 true_time <- rexp(n, rate = true_hazard)

 # Add censoring (approx 30%)
 cens_time <- rexp(n, rate = 0.3)
 obs_time <- pmin(true_time, cens_time)
 obs_status <- as.numeric(true_time <= cens_time)

 # 2. Assign Blocks
 blocks <- sort(rep(1:num_blocks, length.out = n))

 # 3. Induce Linkage Errors (shuffling covariates within blocks)
 linked_x1 <- x1
 linked_x2 <- x2

 for (v in 1:num_blocks) {
  idx <- which(blocks == v)
  n_v <- length(idx)
  mis_size <- round(m_rate * n_v)

  if (mis_size > 0) {
   mis_idx <- sample(idx, mis_size)
   shuff_idx <- sample(mis_idx)
   linked_x1[mis_idx] <- x1[shuff_idx]
   linked_x2[mis_idx] <- x2[shuff_idx]
  }
 }

 dat <- data.frame(
  time = obs_time,
  status = obs_status,
  x1 = linked_x1,
  x2 = linked_x2,
  blocks = blocks
 )

 return(list(data = dat, true_beta = true_beta))
}

context("CoxPH with Exchangeable Linkage Errors (coxphELE)")

test_that("coxphELE engine runs and produces valid output", {
 sim <- generate_ele_cox_data(n = 300, m_rate = 0.1)

 fit <- coxphELE(
  x = as.matrix(sim$data[, c("x1", "x2")]),
  y = sim$data$time,
  cens = 1 - sim$data$status, # Engine uses 1 = censored
  m.rate = 0.1,
  blocks = sim$data$blocks
 )

 expect_s3_class(fit, "coxphELE")
 expect_equal(length(fit$coefficients), 2)
 expect_true(all(!is.na(fit$coefficients)))
 expect_true(all(diag(fit$var) > 0))
})

test_that("S3 Methods (print, summary, vcov, confint) work correctly", {
 sim <- generate_ele_cox_data(n = 200)
 adj <- adjELE(linked.data = sim$data, m.rate = 0.15)
 fit <- plcoxph(Surv(time, status) ~ x1 + x2, adjustment = adj)

 # vcov
 expect_equal(dim(vcov(fit)), c(2, 2))

 # summary
 s_fit <- summary(fit)
 expect_s3_class(s_fit, "summary.coxphELE")
 expect_equal(nrow(s_fit$coefficients), 2)
 expect_equal(colnames(s_fit$coefficients), c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)"))

 # confint
 ci <- confint(fit, level = 0.9)
 expect_equal(dim(ci), c(2, 2))
 expect_equal(colnames(ci), c("5 %", "95 %"))

 # print (capturing output to ensure no errors)
 expect_output(print(fit))
 expect_output(print(s_fit))
})

test_that("Predict method handles lp and risk types", {
 sim <- generate_ele_cox_data(n = 100)
 adj <- adjELE(linked.data = sim$data, m.rate = 0.15)
 fit <- plcoxph(Surv(time, status) ~ x1 + x2, adjustment = adj)

 # Original data
 lp <- predict(fit, type = "lp")
 expect_length(lp, 100)

 risk <- predict(fit, type = "risk")
 expect_true(all(risk > 0))
 expect_equal(risk, exp(lp))

 # New data
 new_dat <- data.frame(x1 = c(0, 1), x2 = c(1, 0))
 lp_new <- predict(fit, newdata = new_dat)
 expect_length(lp_new, 2)
})

test_that("coxphELE handles audit.size for variance inflation", {
 # When audit.size is provided, V1 component of variance is added
 sim <- generate_ele_cox_data(n = 400, m_rate = 0.2, num_blocks = 2)

 # Fit without audit (V1 = 0)
 fit_no_audit <- coxphELE(
  x = as.matrix(sim$data[, c("x1", "x2")]),
  y = sim$data$time,
  cens = 1 - sim$data$status,
  m.rate = 0.2,
  blocks = sim$data$blocks
 )

 # Fit with audit (V1 > 0)
 fit_audit <- coxphELE(
  x = as.matrix(sim$data[, c("x1", "x2")]),
  y = sim$data$time,
  cens = 1 - sim$data$status,
  m.rate = 0.2,
  blocks = sim$data$blocks,
  audit.size = 20 # 20 per block
 )

 # Variances should generally be larger when accounting for m.rate estimation error
 expect_true(all(diag(fit_audit$var) >= diag(fit_no_audit$var)))
})

test_that("fitcoxph.adjELE correctly subsets and aligns data", {
 sim <- generate_ele_cox_data(n = 200)
 # Create some NAs in a variable NOT in the model to test alignment
 sim$data$unused_var <- c(rep(NA, 10), rep(1, 190))

 adj <- adjELE(linked.data = sim$data, m.rate = 0.15)

 # plcoxph internally calls fitcoxph.adjELE
 # Test with a subset to ensure row names/indices align properly
 fit_sub <- plcoxph(Surv(time, status) ~ x1 + x2,
                    data = sim$data,
                    subset = 51:200,
                    adjustment = adj)

 expect_equal(fit_sub$n, 150)
 expect_equal(length(fit_sub$linear.predictors), 150)
})

test_that("Engine throws appropriate errors for bad inputs", {
 sim <- generate_ele_cox_data(n = 50)

 # Mismatch in m.rate and blocks
 expect_error(
  coxphELE(x = as.matrix(sim$data[,3:4]), y = sim$data$time, cens = 1-sim$data$status,
           m.rate = c(0.1, 0.2), blocks = rep(1, 50))
 )

 # Bad response type in wrapper
 adj <- adjELE(linked.data = sim$data, m.rate = 0.1)
 expect_error(
  fitcoxph.adjELE(x = as.matrix(sim$data[,3:4]), y = sim$data$time, adjustment = adj)
 )
})
