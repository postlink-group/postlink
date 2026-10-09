local_edition(3)

# Helper: Create Dummy Data
setup_data <- function() {
 data.frame(
  id = 1:50,
  val = rnorm(50)
 )
}

test_that("print.adjMixBayes: Output for specified data", {
 df <- setup_data()
 adj <- adjMixBayes(linked.data = df)

 # Capture the print output
 out <- capture_output(print(adj))

 # Check for header
 expect_match(out, "Adjustment Object: Bayesian Mixture")

 # Check data summary
 expect_match(out, "Observations:\\s+50")
})

test_that("print.adjMixBayes: Output for NULL data", {
 adj <- adjMixBayes(linked.data = NULL)

 out <- capture_output(print(adj))

 # Check for header
 expect_match(out, "Adjustment Object: Bayesian Mixture")

 # Check status message
 expect_match(out, "Status:\\s+None specified")
})

test_that("print.adjMixBayes: Robustness against corrupted objects", {
 # Simulate an object where the environment exists but is empty
 adj_empty_env <- structure(
  list(data_ref = new.env()),
  class = c("adjMixBayes", "adjustment")
 )

 out <- capture_output(print(adj_empty_env))

 # Should not crash, should report None specified
 expect_match(out, "Status:\\s+None specified")
})

test_that("print.adjMixBayes: reports the match probability model", {
 df <- data.frame(y = rnorm(30), z1 = runif(30),
                  safe = c(rep(TRUE, 10), rep(FALSE, 20)))

 out <- capture_output(print(adjMixBayes(linked.data = df)))
 expect_match(out, "theta ~ beta\\(1,1\\)")
 expect_match(out, "Formula:\\s+~1 \\(constant match probability\\)")
 expect_match(out, "Prior Mismatch Rate:\\s+None specified")
 expect_match(out, "Safe Matches:\\s+None specified")

 adj <- adjMixBayes(linked.data = df, m.rate = 0.2, m.rate.sd = 0.05)
 out <- capture_output(print(adj))
 bp <- postlink:::mrate_to_beta(0.2, 0.05)
 expect_match(out, sprintf("theta ~ beta\\(%s,%s\\) \\[from m.rate\\]",
                           format(bp$alpha, digits = 3), format(bp$beta, digits = 3)))
 expect_match(out, "Prior Mismatch Rate:\\s+0.2 \\(prior SD of theta 0.05\\)")

 adj2 <- adjMixBayes(linked.data = df, m.formula = ~ z1, m.rate = 0.2,
                     safe.matches = safe)
 out <- capture_output(print(adj2))
 expect_match(out, "Formula:\\s+~z1 \\(logistic regression\\)")
 expect_match(out, "gamma intercept ~ normal\\(")
 expect_match(out, "slopes ~ normal\\(0,2.5\\)")
 expect_match(out, "Safe Matches:\\s+10 \\(33.3%\\)")
})
