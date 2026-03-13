#' Methods for Bayesian mixture survival regression fits
#'
#' @description
#' S3 methods for objects returned by \code{survregMixBayes()}, including
#' printing, summarizing fitted models, computing credible intervals and
#' posterior covariance matrices, generating predictions, and pooling Cox
#' regression fits across posterior match classifications.
#'
#' @name mixture_bayessurvreg_methods
#' @keywords internal
NULL

#' Print a survMixBayes model object
#'
#' Prints the model call and posterior mean regression coefficients for the
#' first mixture component of the fitted survival model. In this package,
#' component 1 is interpreted as the correct-match component and component 2
#' as the incorrect-match component.
#'
#' @param x An object of class \code{survMixBayes}.
#' @param digits Minimum number of significant digits to show.
#' @param ... Further arguments (unused).
#'
#' @return The input \code{x}, invisibly.
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
    cat("Coefficients (posterior means, component 1 = correct-match):\n")
    print(round(colMeans(b), digits = digits))
  }
  invisible(x)
}

#' Summary method for survMixBayes models
#'
#' Computes posterior summaries for the regression coefficients, mixing weight,
#' and component-specific distribution parameters in a fitted
#' \code{survMixBayes} model. Throughout, component 1 is interpreted as the
#' correct-match component and component 2 as the incorrect-match component.
#'
#' @param object An object of class \code{survMixBayes}.
#' @param probs Numeric vector of probabilities used to compute posterior
#'   quantiles for the model parameters. The default,
#'   \code{c(0.025, 0.5, 0.975)}, gives a posterior median and a 95\%
#'   credible interval.
#' @param ... Further arguments (unused).
#'
#' @return An object of class \code{summary.survMixBayes} containing posterior
#'   quantile summaries for the regression coefficients in both mixture
#'   components, the mixing weight, and any family-specific distribution
#'   parameters included in the fitted model. Component 1 corresponds to the
#'   correct-match component and component 2 to the incorrect-match component.
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
#' # Generate and print posterior summaries
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

  fmt(x$coef1, "Coefficients (component 1 = correct-match)")
  fmt(x$coef2, "Coefficients (component 2 = incorrect-match)")

  cat("\nTheta (mix weight for component 1 = correct-match):\n")
  print(round(x$theta, digits = digits))

  fmt(x$shape1, "Shape (component 1 = correct-match)")
  fmt(x$shape2, "Shape (component 2 = incorrect-match)")
  fmt(x$scale1, "Scale (component 1 = correct-match)")
  fmt(x$scale2, "Scale (component 2 = incorrect-match)")

  invisible(x)
}

#' Credible intervals for parameters from a survMixBayes fit
#'
#' Computes posterior credible intervals for the regression coefficients,
#' mixing weight, and family-specific distribution parameters from a fitted
#' \code{survMixBayes} model. The returned intervals are organized by mixture
#' component and parameter type. Component 1 corresponds to the correct-match
#' component and component 2 to the incorrect-match component.
#'
#' @param object An object of class \code{survMixBayes}.
#' @param parm Optional character vector selecting which interval blocks to
#'   return. For example, \code{"theta"} returns only the credible interval for
#'   the mixing weight. If \code{NULL}, credible intervals are returned for all
#'   available parameter blocks.
#' @param level Probability level for the credible intervals. Defaults to
#'   \code{0.95}.
#' @param ... Further arguments (unused).
#'
#' @return A named list of credible intervals. Elements \code{coef1} and
#'   \code{coef2} are matrices with one row per regression coefficient and two
#'   columns giving the lower and upper interval bounds for components 1 and 2,
#'   respectively, where component 1 is the correct-match component and
#'   component 2 is the incorrect-match component. Elements such as \code{theta}, \code{shape1},
#'   \code{shape2}, \code{scale1}, and \code{scale2} are numeric vectors of
#'   length 2 containing the lower and upper credible interval bounds for the
#'   corresponding scalar parameters.
#'
#' @examples
#' \donttest{
#' # Simulate simple Weibull survival data
#' set.seed(403)
#' n <- 100
#'
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' X <- cbind(x1 = x1, x2 = x2)
#'
#' # Observed survival response: time and event indicator
#' y <- cbind(
#'   time = rweibull(n, shape = 1.2, scale = exp(x1)),
#'   event = rbinom(n, 1, 0.8)
#' )
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

#' Posterior covariance matrix for survMixBayes coefficients
#'
#' Returns the empirical posterior covariance matrix of the regression
#' coefficients for component 1 of a fitted \code{survMixBayes} model.
#' In this package, component 1 is interpreted as the correct-match component.
#'
#' @param object A \code{survMixBayes} model object.
#' @param ... Further arguments (unused).
#'
#' @return Posterior covariance matrix of the regression coefficients for
#'   component 1, interpreted as the correct-match component.
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
#' Returns posterior mean linear predictors for each component. Component 1 is
#' interpreted as the correct-match component and component 2 as the incorrect-match
#' component.
#'
#' @param object An object of class \code{survMixBayes}.
#' @param newdata Optional design matrix for prediction. If \code{NULL}, uses \code{object$X} if present.
#' @param ... Additional arguments (unused).
#'
#' @return A list with posterior mean linear predictors for each component:
#'   \code{component1} for the correct-match component and \code{component2}
#'   for the incorrect-match component.
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

#' Pool regression fits across posterior draws of correct-match classifications
#'
#' @description
#' Use posterior draws of the latent match indicators from \code{survregMixBayes()}
#' to repeatedly identify which records are treated as correct matches, refit a
#' Cox proportional hazards model on those records, and pool the resulting
#' estimates using multiple-imputation pooling rules.
#'
#' Each retained posterior draw defines one subset of records classified as
#' correct matches. The function fits the specified \code{survival::coxph()}
#' model to that subset, extracts the estimated coefficients and covariance
#' matrix, and combines the results across draws using Rubin's rules.
#'
#' @param object A \code{survMixBayes} model object containing posterior draws of
#'   the latent match indicators.
#' @param data A data.frame with all candidate records in the same row order as used in the model.
#' @param formula Model formula for refitting on each draw (required), typically
#'   of the form \code{survival::Surv(time, event) ~ ...}.
#' @param min_n Minimum number of records required to fit the model for a given
#'   posterior draw. The default is \code{p + 2}, where \code{p} is the number
#'   of non-intercept columns in the model matrix.
#' @param quietly If \code{TRUE}, draws that lead to fitting errors are skipped
#'   without printing the full error message.
#' @param ties Method for handling tied event times in \code{survival::coxph()}.
#'   Default is \code{"efron"}.
#' @param ... Additional arguments passed to \code{survival::coxph()}.
#'
#' @return An object of class \code{c("mi_link_pool_survreg", "mi_link_pool")}
#'   containing pooled coefficient estimates, standard errors, confidence
#'   intervals, and related summary information.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' set.seed(406)
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' X <- cbind(x1, x2)
#' y <- cbind(
#'   time = rweibull(n, shape = 1.2, scale = exp(x1)),
#'   event = rbinom(n, 1, 0.8)
#' )
#'
#' fit <- survregMixBayes(
#'   X = X, y = y, dist = "weibull",
#'   control = list(iterations = 100, burnin.iterations = 50, seed = 406)
#' )
#'
#' dat <- data.frame(time = y[, "time"], event = y[, "event"], x1 = x1, x2 = x2)
#'
#' pooled_obj <- mi_with(
#'   object = fit,
#'   data = dat,
#'   formula = survival::Surv(time, event) ~ x1 + x2
#' )
#'
#' print(pooled_obj)
#' }
#'
#' @export
#' @method mi_with survMixBayes
#' @importFrom stats coef vcov cov qt model.frame model.matrix
#' @importFrom survival coxph Surv
mi_with.survMixBayes <- function(object, data, formula,
                                 min_n = NULL, quietly = TRUE,
                                 ties = "efron", ...) {

 if (missing(formula) || is.null(formula)) {
  if (!is.null(object$call$formula)) {
   formula <- eval(object$call$formula, envir = parent.frame())
  } else {
   stop(
    "Formula not found: please provide it explicitly or ensure object$call$formula exists.",
    call. = FALSE
   )
  }
 }

 if (!inherits(formula, "formula")) {
  stop("`formula` must be a valid formula object.", call. = FALSE)
 }

 if (missing(data) || is.null(data) || !is.data.frame(data)) {
  stop("`data` must be provided as a data.frame.", call. = FALSE)
 }

 z_samples <- object$m_samples

 collapse_z_g <- function(z_samples, g) {
  if (is.list(z_samples)) {
   parts <- lapply(z_samples, function(Z) {
    if (is.null(dim(Z))) {
     stop("Each list element of `z_samples` must be a matrix-like object.", call. = FALSE)
    }
    Z[g, ]
   })
   as.integer(do.call(c, parts))
  } else {
   if (is.null(dim(z_samples)) || length(dim(z_samples)) != 2L) {
    stop("`z_samples` must be an S x N matrix, or a list of such matrices.", call. = FALSE)
   }
   as.integer(z_samples[g, ])
  }
 }

 fit_once <- function(df) {
  fit <- survival::coxph(formula = formula, data = df, ties = ties, ...)
  list(coef = stats::coef(fit), vcov = stats::vcov(fit))
 }

 mf <- stats::model.frame(formula, data = data)
 Xtmp <- stats::model.matrix(formula, mf)
 if ("(Intercept)" %in% colnames(Xtmp)) {
  Xtmp <- Xtmp[, colnames(Xtmp) != "(Intercept)", drop = FALSE]
 }
 p <- ncol(Xtmp)
 if (is.null(min_n)) min_n <- p + 2L

 S <- if (is.list(z_samples)) nrow(z_samples[[1]]) else nrow(z_samples)
 if (length(S) == 0L || is.null(S) || !is.finite(S)) {
  stop("Unable to infer the number of posterior draws (S).", call. = FALSE)
 }

 coefs_list <- list()
 vcovs_list <- list()
 kept <- logical(S)

 for (s in seq_len(S)) {
  zs <- collapse_z_g(z_samples, s)

  if (!all(zs %in% c(1L, 2L))) {
   stop("`m_samples` must contain only 1/2 indicators.", call. = FALSE)
  }

  idx <- which(zs == 1L)
  if (length(idx) < min_n) next

  dat_s <- data[idx, , drop = FALSE]

  res <- try(fit_once(dat_s), silent = quietly)
  if (!inherits(res, "try-error") &&
      is.numeric(res$coef) &&
      is.matrix(res$vcov) &&
      length(res$coef) > 0L) {
   coefs_list[[length(coefs_list) + 1L]] <- res$coef
   vcovs_list[[length(vcovs_list) + 1L]] <- res$vcov
   kept[s] <- TRUE
  }
 }

 m <- length(coefs_list)
 if (m == 0L) {
  stop(
   "No valid imputations: all draws failed or had too few component-1 records.",
   call. = FALSE
  )
 }

 parnames <- names(coefs_list[[1L]])
 coefs_mat <- do.call(rbind, lapply(coefs_list, function(b) b[parnames]))

 Ubar <- Reduce("+", vcovs_list) / m
 B <- stats::cov(coefs_mat)
 Tmat <- Ubar + (1 + 1 / m) * B

 qbar <- colMeans(coefs_mat)
 se <- sqrt(diag(Tmat))
 lambda <- diag((1 + 1 / m) * B) / diag(Tmat)

 r <- diag((1 + 1 / m) * B) / diag(Ubar)
 r[!is.finite(r)] <- 0
 dfold <- (m - 1) * (1 + 1 / pmax(r, .Machine$double.eps))^2
 df <- pmax(dfold, 3)

 tcrit <- stats::qt(0.975, df = df)
 lwr <- qbar - tcrit * se
 upr <- qbar + tcrit * se
 ci95 <- cbind(lwr = lwr, upr = upr)
 rownames(ci95) <- names(qbar)

 out <- list(
  m          = m,
  coef       = qbar,
  vcov       = Tmat,
  se         = se,
  ci95       = ci95,
  Ubar       = Ubar,
  B          = B,
  lambda     = lambda,
  df         = df,
  kept_draws = which(kept),
  dist       = object$dist,
  refit      = "coxph",
  call       = object$call
 )

 class(out) <- c("mi_link_pool_survreg", "mi_link_pool")
 out
}

#' Print pooled Cox regression results
#'
#' @param x An object of class \code{mi_link_pool_survreg}, typically returned by
#'   \code{mi_with()} for a \code{survMixBayes} fit.
#' @param digits the number of significant digits to print.
#' @param ... further arguments (unused).
#'
#' @return The input \code{x}, invisibly.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' set.seed(407)
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' X <- cbind(x1, x2)
#' y <- cbind(
#'   time = rweibull(n, shape = 1.2, scale = exp(x1)),
#'   event = rbinom(n, 1, 0.8)
#' )
#'
#' fit <- survregMixBayes(
#'   X = X, y = y, dist = "weibull",
#'   control = list(iterations = 100, burnin.iterations = 50, seed = 407)
#' )
#'
#' dat <- data.frame(time = y[, "time"], event = y[, "event"], x1 = x1, x2 = x2)
#'
#' pooled_obj <- mi_with(
#'   object = fit,
#'   data = dat,
#'   formula = survival::Surv(time, event) ~ x1 + x2
#' )
#'
#' print(pooled_obj, digits = 4)
#' }
#'
#' @export
#' @method print mi_link_pool_survreg
print.mi_link_pool_survreg <- function(x,
                                       digits = max(3L, getOption("digits") - 2L),
                                       ...) {
 cat("Pooled Cox regression results across posterior match classifications:\n")
 cat("  Retained imputations (m):", x$m, "\n")
 cat("  Mixture model distribution:", x$dist, "\n")
 cat("  Refit model: coxph\n\n")

 tab <- cbind(
  Estimate  = x$coef,
  Std.Error = x$se,
  CI.lwr    = x$ci95[, "lwr"],
  CI.upr    = x$ci95[, "upr"],
  df        = x$df
 )

 print(round(tab, digits))
 invisible(x)
}
