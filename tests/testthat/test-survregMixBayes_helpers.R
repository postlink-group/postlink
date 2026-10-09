# tests/testthat/test-survregMixBayes_helpers.R
# Unit tests for internal Bayesian survival mixture helper utilities
# These tests are fast and do not run MCMC.

local_edition(3)

test_that(".validate_survreg_dist lowercases valid input and rejects invalid shape", {
 expect_equal(postlink:::.validate_survreg_dist("weibull"), "weibull")
 expect_equal(postlink:::.validate_survreg_dist("Gamma"), "gamma")
 expect_equal(postlink:::.validate_survreg_dist("WEIBULL"), "weibull")

 expect_error(
  postlink:::.validate_survreg_dist(),
  "`dist` must be a single character string."
 )

 expect_error(
  postlink:::.validate_survreg_dist(1),
  "`dist` must be a single character string."
 )

 expect_error(
  postlink:::.validate_survreg_dist(c("gamma", "weibull")),
  "`dist` must be a single character string."
 )
})

test_that(".normalize_surv_y works for list input", {
 y <- list(
  time = c(1.2, 2.5, 3.1),
  event = c(1, 0, 1)
 )

 out <- postlink:::.normalize_surv_y(y)

 expect_type(out, "list")
 expect_equal(names(out), c("time", "event"))
 expect_equal(out$time, c(1.2, 2.5, 3.1))
 expect_equal(out$event, c(1L, 0L, 1L))

 # event codes other than 0/1 are refused instead of being turned into events
 expect_error(postlink:::.normalize_surv_y(list(time = c(1.2, 2.5), event = c(1, 3))), "coded 0")
 # logical events are accepted
 expect_equal(postlink:::.normalize_surv_y(list(time = c(1, 2), event = c(TRUE, FALSE)))$event, c(1L, 0L))
})

test_that(".normalize_surv_y works for matrix input", {
 y <- cbind(
  time = c(0.5, 1.5, 2.5),
  event = c(0, 1, 1)
 )

 out <- postlink:::.normalize_surv_y(y)

 expect_equal(out$time, c(0.5, 1.5, 2.5))
 expect_equal(out$event, c(0L, 1L, 1L))

 expect_error(postlink:::.normalize_surv_y(cbind(c(0.5, 1.5), c(1, 2))), "coded 0")
 expect_error(postlink:::.normalize_surv_y(cbind(c(0.5, 1.5), c(1, 0), c(1, 1))), "2-column")
 expect_equal(postlink:::.normalize_surv_y(survival::Surv(c(0.5, 1.5), c(1, 0)))$event, c(1L, 0L))
 expect_error(postlink:::.normalize_surv_y(survival::Surv(c(0.5, 1.5), c(1, 0), type = "left")),
              "right-censored")
})

test_that(".normalize_surv_y rejects invalid inputs", {
 expect_error(
  postlink:::.normalize_surv_y(c(1, 2, 3)),
  "`y` must be a 2-column matrix"
 )

 expect_error(
  postlink:::.normalize_surv_y(matrix(c(1, 0, 2), ncol = 1)),
  "`y` must be a 2-column matrix"
 )

 expect_error(
  postlink:::.normalize_surv_y(list(time = c(1, 2, 3))),
  "`y` must be a 2-column matrix"
 )

 expect_error(
  postlink:::.normalize_surv_y(list(time = c(1, 0, 2), event = c(1, 0, 1))),
  "Survival times must be positive and finite."
 )

 expect_error(
  postlink:::.normalize_surv_y(list(time = c(1, -2, 3), event = c(1, 0, 1))),
  "Survival times must be positive and finite."
 )

 expect_error(
  postlink:::.normalize_surv_y(list(time = c(1, Inf, 3), event = c(1, 0, 1))),
  "Survival times must be positive and finite."
 )
})
