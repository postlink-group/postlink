# Methods for survMixBayes objects
# S3 methods for the Bayesian two-component survival mixture model.

#' Print a survMixBayes object
#'
#' Prints the model call and posterior mean regression coefficients for component 1.
#'
#' @param x An object of class \code{survMixBayes}.
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments (unused).
#'
#' @return The object \code{x} (invisibly).
#'
#' @examples
#' \donttest{
#' # 1. Simulate fast toy data
#' set.seed(401)
#' n <- 100
#' X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
#' y <- cbind(time = rweibull(n, shape = 1.5, scale = exp(0.5 * X[, 1])),
#'            event = rbinom(n, 1, 0.8))
#'
#' # 2. Fit the model with artificially low iterations for speed
#' fit <- survregMixBayes(
#'   X = X, y = y, dist = "weibull",
#'   control = list(iterations = 100, burnin.iterations = 50, seed = 401)
#' )
#'
#' # 3. Print the model object
#' print(fit)
#' }
#'
#' @export
#' @method print survMixBayes
print.survMixBayes <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Bayesian two-component mixture survival regression\n")
  cat("Family:", x$dist, "\n\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  b <- x$estimates$coefficients
  if (is.matrix(b)) {
    cat("Coefficients (posterior means, component 1):\n")
    print(round(colMeans(b), digits = digits))
  }
  invisible(x)
}

#' Summarize a survMixBayes object
#'
#' Computes posterior summaries for coefficients and key distribution parameters.
#'
#' @param object An object of class \code{survMixBayes}.
#' @param probs Numeric vector of quantile probabilities.
#' @param ... Additional arguments (unused).
#'
#' @return An object of class \code{summary.survMixBayes}.
#'
#' @examples
#' \donttest{
#' # Simulate data
#' set.seed(402)
#' n <- 100
#' X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
#' y <- cbind(time = rweibull(n, shape = 1.2, scale = exp(X[, 1])),
#'            event = rbinom(n, 1, 0.8))
#'
#' # Fit the model
#' fit <- survregMixBayes(
#'   X = X, y = y, dist = "weibull",
#'   control = list(iterations = 100, burnin.iterations = 50, seed = 402)
#' )
#'
#' # Generate and print the comprehensive summary
#' fit_summary <- summary(fit, probs = c(0.025, 0.5, 0.975))
#' print(fit_summary)
#' }
#'
#' @export
#' @method summary survMixBayes
summary.survMixBayes <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  qfun <- function(x) stats::quantile(x, probs = probs, na.rm = TRUE)
  b1 <- object$estimates$coefficients
  b2 <- object$estimates$m.coefficients
  if (!is.matrix(b1) || !is.matrix(b2)) stop("Expected coefficient draws as matrices.", call. = FALSE)

  s <- list(
    call = object$call,
    family = object$dist,
    coef1 = apply(b1, 2, qfun),
    coef2 = apply(b2, 2, qfun),
    theta = qfun(object$estimates$theta)
  )

  if (!is.null(object$estimates$shape))  s$shape1 <- qfun(object$estimates$shape)
  if (!is.null(object$estimates$m.shape)) s$shape2 <- qfun(object$estimates$m.shape)
  if (!is.null(object$estimates$scale))  s$scale1 <- qfun(object$estimates$scale)
  if (!is.null(object$estimates$m.scale)) s$scale2 <- qfun(object$estimates$m.scale)

  class(s) <- "summary.survMixBayes"
  s
}

#' @noRd
#' @export
print.summary.survMixBayes <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Summary of Bayesian mixture survival regression\n")
  if (!is.null(x$call)) {
    cat("\nCall:\n")
    print(x$call)
  }
  cat("\nFamily:", x$family, "\n")

  fmt <- function(mat, nm) {
    if (is.null(mat)) return(invisible(NULL))
    cat("\n", nm, " (quantiles):\n", sep = "")
    print(round(t(mat), digits = digits))
  }

  fmt(x$coef1, "Coefficients (component 1)")
  fmt(x$coef2, "Coefficients (component 2)")

  cat("\nTheta (mix weight for component 1):\n")
  print(round(x$theta, digits = digits))

  fmt(x$shape1, "Shape (component 1)")
  fmt(x$shape2, "Shape (component 2)")
  fmt(x$scale1, "Scale (component 1)")
  fmt(x$scale2, "Scale (component 2)")

  invisible(x)
}

#' Confidence intervals for survMixBayes parameters
#'
#' Returns posterior credible intervals for regression coefficients and key parameters.
#'
#' @param object An object of class \code{survMixBayes}.
#' @param parm Optional. Character vector selecting which elements of the returned
#'   list to keep. If \code{NULL}, all intervals are returned.
#' @param level Confidence level.
#' @param ... Additional arguments (unused).
#' @return A named list of credible intervals. Elements `coef1` and `coef2`
#'   are `p x 2` matrices (lower/upper) for component-specific regression
#'   coefficients when available. Elements `theta`, `shape1`, `shape2`,
#'   `scale1`, `scale2` are numeric vectors of length 2 giving lower/upper
#'   credible intervals for scalar parameters.
#'
#' @examples
#' \donttest{
#' # Simulate data
#' set.seed(403)
#' n <- 100
#' X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
#' y <- cbind(time = rweibull(n, shape = 1.2, scale = exp(X[, 1])),
#'            event = rbinom(n, 1, 0.8))
#'
#' # Fit the model
#' fit <- survregMixBayes(
#'   X = X, y = y, dist = "weibull",
#'   control = list(iterations = 100, burnin.iterations = 50, seed = 403)
#' )
#'
#' # Calculate 95% credible intervals for all parameters
#' confint(fit, level = 0.95)
#'
#' # Extract credible intervals specifically for the mixing weight
#' confint(fit, parm = "theta", level = 0.90)
#' }
#'
#' @export
#' @method confint survMixBayes
confint.survMixBayes <- function(object, parm = NULL, level = 0.95, ...) {
 alpha <- (1 - level) / 2
 probs <- c(alpha, 1 - alpha)

 out <- list()

 b1 <- object$estimates$coefficients
 b2 <- object$estimates$m.coefficients
 if (is.matrix(b1)) out$coef1 <- t(apply(b1, 2, stats::quantile, probs = probs, na.rm = TRUE))
 if (is.matrix(b2)) out$coef2 <- t(apply(b2, 2, stats::quantile, probs = probs, na.rm = TRUE))

 out$theta <- stats::quantile(object$estimates$theta, probs = probs, na.rm = TRUE)

 if (!is.null(object$estimates$shape))   out$shape1 <- stats::quantile(object$estimates$shape,   probs = probs, na.rm = TRUE)
 if (!is.null(object$estimates$m.shape)) out$shape2 <- stats::quantile(object$estimates$m.shape, probs = probs, na.rm = TRUE)
 if (!is.null(object$estimates$scale))   out$scale1 <- stats::quantile(object$estimates$scale,   probs = probs, na.rm = TRUE)
 if (!is.null(object$estimates$m.scale)) out$scale2 <- stats::quantile(object$estimates$m.scale, probs = probs, na.rm = TRUE)

 # Optional filtering (simple and safe): allow selecting list elements by name
 if (!is.null(parm)) {
  if (is.character(parm)) out <- out[names(out) %in% parm]
 }

 out
}

#' Variance-covariance for survMixBayes coefficients
#'
#' Returns an empirical posterior covariance matrix of component 1 regression coefficients.
#'
#' @param object An object of class \code{survMixBayes}.
#' @param ... Additional arguments (unused).
#'
#' @return A covariance matrix.
#'
#' @examples
#' \donttest{
#' # Simulate data
#' set.seed(404)
#' n <- 100
#' X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
#' y <- cbind(time = rweibull(n, shape = 1.2, scale = exp(X[, 1])),
#'            event = rbinom(n, 1, 0.8))
#'
#' # Fit the model
#' fit <- survregMixBayes(
#'   X = X, y = y, dist = "weibull",
#'   control = list(iterations = 100, burnin.iterations = 50, seed = 404)
#' )
#'
#' # Extract the empirical posterior covariance matrix for component 1
#' vcov_mat <- vcov(fit)
#' print(vcov_mat)
#' }
#'
#' @export
#' @method vcov survMixBayes
vcov.survMixBayes <- function(object, ...) {
  b <- object$estimates$coefficients
  if (!is.matrix(b)) stop("Coefficient draws are not available.", call. = FALSE)
  stats::cov(b)
}

#' Predict for survMixBayes
#'
#' Returns posterior mean linear predictors for each component.
#'
#' @param object An object of class \code{survMixBayes}.
#' @param newdata Optional design matrix for prediction. If \code{NULL}, uses \code{object$X} if present.
#' @param ... Additional arguments (unused).
#'
#' @return A list with posterior mean linear predictors for each component.
#'
#' @examples
#' \donttest{
#' # Simulate data
#' set.seed(405)
#' n <- 100
#' X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
#' y <- cbind(time = rweibull(n, shape = 1.2, scale = exp(X[, 1])),
#'            event = rbinom(n, 1, 0.8))
#'
#' # Fit the model
#' fit <- survregMixBayes(
#'   X = X, y = y, dist = "weibull",
#'   control = list(iterations = 100, burnin.iterations = 50, seed = 405)
#' )
#'
#' # Create a new design matrix for prediction
#' X_new <- matrix(c(0, 1, -1, 0.5), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
#'
#' # Predict posterior mean linear predictors for each latent component
#' preds <- predict(fit, newdata = X_new)
#' print(preds$component1)
#' print(preds$component2)
#' }
#'
#' @export
#' @method predict survMixBayes
predict.survMixBayes <- function(object, newdata = NULL, ...) {
  X <- newdata
  if (is.null(X)) {
    if (!is.null(object$X)) X <- object$X
  }
  if (is.null(X)) stop("`newdata` must be provided (object does not store X).", call. = FALSE)
  if (!is.matrix(X) || !is.numeric(X)) stop("`newdata` must be a numeric matrix.", call. = FALSE)

  b1 <- object$estimates$coefficients
  b2 <- object$estimates$m.coefficients
  mu1 <- as.numeric(X %*% colMeans(b1))
  mu2 <- as.numeric(X %*% colMeans(b2))
  list(component1 = mu1, component2 = mu2)
}

#' Pooling method for survMixBayes objects
#'
#' Performs posterior allocation based pooling for a \code{survMixBayes} fit.
#'
#' @param object An object of class \code{survMixBayes}.
#' @param ... Additional arguments (unused).
#'
#' @return A pooled object of class \code{mi_link_pool_survreg}.
#'
#' @examples
#' \donttest{
#' # Simulate data
#' set.seed(406)
#' n <- 100
#' X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
#' y <- cbind(time = rweibull(n, shape = 1.2, scale = exp(X[, 1])),
#'            event = rbinom(n, 1, 0.8))
#'
#' # Fit the model
#' fit <- survregMixBayes(
#'   X = X, y = y, dist = "weibull",
#'   control = list(iterations = 100, burnin.iterations = 50, seed = 406)
#' )
#'
#' # Perform posterior allocation pooling
#' pooled_obj <- mi_with(fit)
#'
#' # View the distribution of posterior probabilities of belonging to component 1
#' summary(pooled_obj$p_component1)
#' }
#'
#' @export
#' @method mi_with survMixBayes
mi_with.survMixBayes <- function(object, ...) {

 z <- object$m_samples
 if (!is.matrix(z)) stop("Expected `m_samples` as a matrix (S x N).")

 # Posterior probability of being in component 1
 p1 <- colMeans(z == 1L, na.rm = TRUE)

 pooled <- list(
  call = object$call,
  dist = object$dist,
  p_component1 = p1,
  estimates = object$estimates
 )
 class(pooled) <- "mi_link_pool_survreg"
 pooled
}

#' Print method for pooled survreg MI object
#'
#' @param x An object of class \code{mi_link_pool_survreg}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' # Simulate data
#' set.seed(407)
#' n <- 100
#' X <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
#' y <- cbind(time = rweibull(n, shape = 1.2, scale = exp(X[, 1])),
#'            event = rbinom(n, 1, 0.8))
#'
#' # Fit the model
#' fit <- survregMixBayes(
#'   X = X, y = y, dist = "weibull",
#'   control = list(iterations = 100, burnin.iterations = 50, seed = 407)
#' )
#'
#' # Create the pooled object via mi_with
#' pooled_obj <- mi_with(fit)
#'
#' # Explicitly call the print method for the pooled object
#' print(pooled_obj)
#' }
#'
#' @export
#' @method print mi_link_pool_survreg
print.mi_link_pool_survreg <- function(x, ...) {

 cat("Pooled posterior allocation (component 1):\n")
 print(summary(x$p_component1))
 invisible(x)
}
