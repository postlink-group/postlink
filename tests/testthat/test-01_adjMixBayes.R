local_edition(3)

# Helper: Create Dummy Data
setup_data <- function() {
 data.frame(
  id = 1:50,
  val = rnorm(50)
 )
}

test_that("adjMixBayes: Basic valid construction works", {
 df <- setup_data()

 # Initialize with dataframe
 adj <- adjMixBayes(linked.data = df)

 # Check Class
 expect_s3_class(adj, "adjMixBayes")
 expect_s3_class(adj, "adjustment")

 # Check Structure
 expect_true(is.environment(adj$data_ref))
 expect_identical(adj$data_ref$data, df)
})

test_that("adjMixBayes: Reference semantics (lightweight container) are preserved", {
 df <- setup_data()
 adj <- adjMixBayes(linked.data = df)

 # Ensure data is stored in an environment, not directly in the list
 expect_type(adj$data_ref, "environment")

 # Modify the environment data manually to prove it's a reference
 # (In a real scenario, this helps avoid copying the dataframe)
 adj$data_ref$data$new_col <- 1
 expect_true("new_col" %in% names(adj$data_ref$data))
})

test_that("adjMixBayes: Input validation enforces types", {
 # Error on character string
 expect_error(adjMixBayes(linked.data = "invalid_string"),
              "must be a data.frame, list, or environment")

 # Error on numeric vector
 expect_error(adjMixBayes(linked.data = c(1, 2, 3)),
              "must be a data.frame, list, or environment")
})

test_that("adjMixBayes: Handling of coercible types (List)", {
 df <- setup_data()
 df_list <- as.list(df)

 # Should successfully coerce list to data.frame
 adj <- adjMixBayes(linked.data = df_list)

 expect_s3_class(adj$data_ref$data, "data.frame")
 expect_equal(nrow(adj$data_ref$data), 50)
})

test_that("adjMixBayes: Handling of NULL input", {
 # Allow NULL construction
 adj <- adjMixBayes(linked.data = NULL)

 expect_s3_class(adj, "adjMixBayes")
 expect_true(is.environment(adj$data_ref))
 expect_null(adj$data_ref$data)
})

# --- Linkage information: m.formula, m.rate, m.rate.sd, safe.matches ---

test_that("adjMixBayes: m.formula validation", {
 df <- data.frame(y = rnorm(20), z1 = rnorm(20), z2 = rnorm(20))

 # Valid one-sided formula
 adj <- adjMixBayes(linked.data = df, m.formula = ~ z1 + z2)
 expect_equal(deparse(adj$m.formula), "~z1 + z2")

 # Default is ~1 (intercept only); NULL is treated the same way
 adj2 <- adjMixBayes(linked.data = df)
 expect_equal(deparse(adj2$m.formula), "~1")
 expect_true(postlink:::.is_intercept_only_formula(adj2$m.formula))
 expect_equal(deparse(adjMixBayes(linked.data = df, m.formula = NULL)$m.formula), "~1")

 # Error: two-sided formula
 expect_error(adjMixBayes(linked.data = df, m.formula = y ~ z1),
              "one-sided formula")

 # Error: not a formula
 expect_error(adjMixBayes(linked.data = df, m.formula = "~ z1"),
              "must be a formula")

 # Error: variable not in data; '.' unsupported
 expect_error(adjMixBayes(linked.data = df, m.formula = ~ missing_var),
              "not found in 'linked.data'")
 expect_error(adjMixBayes(linked.data = df, m.formula = ~ .), "not supported")
})

test_that("adjMixBayes: m.rate and m.rate.sd validation", {
 df <- setup_data()

 adj <- adjMixBayes(linked.data = df, m.rate = 0.3)
 expect_equal(adj$m.rate, 0.3)
 expect_equal(adj$m.rate.sd, 0.1)

 adj2 <- adjMixBayes(linked.data = df)
 expect_null(adj2$m.rate)

 expect_error(adjMixBayes(linked.data = df, m.rate = 0), "strictly between 0 and 1")
 expect_error(adjMixBayes(linked.data = df, m.rate = 1), "strictly between 0 and 1")
 expect_error(adjMixBayes(linked.data = df, m.rate = -0.1), "strictly between 0 and 1")
 expect_error(adjMixBayes(linked.data = df, m.rate = "0.3"), "single numeric value")

 # probability-scale SD must be positive and compatible with a Beta prior
 expect_error(adjMixBayes(linked.data = df, m.rate = 0.3, m.rate.sd = 0), "positive")
 expect_error(adjMixBayes(linked.data = df, m.rate = 0.3, m.rate.sd = NULL), "positive")
 expect_error(adjMixBayes(linked.data = df, m.rate = 0.05, m.rate.sd = 0.3),
              "strictly less than")
 adj3 <- adjMixBayes(linked.data = df, m.rate = 0.05, m.rate.sd = 0.02)
 expect_equal(adj3$m.rate.sd, 0.02)
})

test_that("adjMixBayes: safe.matches validation and NSE resolution", {
 df <- data.frame(y = rnorm(20), safe = c(rep(TRUE, 5), rep(FALSE, 15)),
                  safe01 = c(rep(1, 5), rep(0, 15)))

 # Direct logical vector
 adj <- adjMixBayes(linked.data = df, safe.matches = df$safe)
 expect_equal(sum(adj$safe.matches), 5)
 expect_length(adj$safe.matches, 20)
 expect_type(adj$safe.matches, "logical")

 # NSE resolution from linked.data, logical or 0/1
 adj2 <- adjMixBayes(linked.data = df, safe.matches = safe)
 expect_equal(sum(adj2$safe.matches), 5)
 adj01 <- adjMixBayes(linked.data = df, safe.matches = safe01)
 expect_identical(adj01$safe.matches, df$safe)

 # Standard evaluation from the calling environment
 my_safe <- df$safe
 adj3 <- adjMixBayes(linked.data = df, safe.matches = my_safe)
 expect_equal(sum(adj3$safe.matches), 5)

 # NULL is fine (default)
 expect_null(adjMixBayes(linked.data = df)$safe.matches)

 # Errors: unknown object, wrong type, NA, wrong length
 expect_error(adjMixBayes(linked.data = df, safe.matches = not_a_column),
              "Could not find object")
 expect_error(adjMixBayes(linked.data = df, safe.matches = y), "logical \\(or 0/1\\) vector")
 expect_error(adjMixBayes(linked.data = df, safe.matches = c(TRUE, NA, rep(FALSE, 18))),
              "NA values")
 expect_error(adjMixBayes(linked.data = df, safe.matches = c(TRUE, FALSE)),
              "same length")
})

test_that("adjMixBayes: combined parameters stored correctly", {
 df <- data.frame(y = rnorm(30), z1 = runif(30),
                  safe = c(rep(TRUE, 10), rep(FALSE, 20)))

 adj <- adjMixBayes(
  linked.data = df,
  priors = list(theta = "beta(70, 30)"),
  m.formula = ~ z1,
  m.rate = 0.2,
  safe.matches = safe
 )

 expect_s3_class(adj, "adjMixBayes")
 expect_equal(deparse(adj$m.formula), "~z1")
 expect_equal(adj$m.rate, 0.2)
 expect_equal(sum(adj$safe.matches), 10)
 expect_equal(adj$priors$theta, "beta(70, 30)")
 expect_error(adjMixBayes(linked.data = df, priors = list("beta(1,1)")), "named list")
})
