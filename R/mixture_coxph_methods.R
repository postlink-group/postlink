#' Extract Variance-Covariance Matrix from a coxphMixture Object
#'
#' @description
#' Extracts the variance-covariance matrix of the main parameters from a fitted 
#' \code{coxphMixture} object. The matrix is estimated using Louis' method (1982) 
#' to account for the missing data structure (latent match status) inherent in the 
#' mixture model.
#'
#' @param object An object of class \code{coxphMixture}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A matrix of the estimated covariances between the parameter estimates. 
#' The rows and columns correspond to the outcome model coefficients (\code{beta}) 
#' and the mismatch model coefficients (\code{gamma}).
#'
#' @references
#' Louis, T. A. (1982). Finding the observed information matrix when using the 
#' EM algorithm. \emph{Journal of the Royal Statistical Society: Series B 
#' (Methodological)}, 44(2), 226-233.
#'
#' @export
vcov.coxphMixture <- function(object, ...) {
 return(object$var)
}

#' Confidence Intervals for coxphMixture Objects
#'
#' @description
#' Computes Wald confidence intervals for the coefficients of the outcome model 
#' and the mismatch indicator model.
#'
#' @param object An object of class \code{coxphMixture}.
#' @param parm A specification of which parameters are to be given confidence intervals, 
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level The confidence level required (default is 0.95).
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' The intervals are calculated based on the variance-covariance matrix returned 
#' by \code{\link{vcov.coxphMixture}}, using the standard normal approximation: 
#' \code{Estimate +/- z_crit * SE}.
#'
#' @return A matrix with columns giving lower and upper confidence limits for each parameter.
#'
#' @export
confint.coxphMixture <- function(object, parm, level = 0.95, ...) {
 # Aggregate all coefficients
 coefs <- c(object$coefficients, object$m.coefficients)
 
 # Ensure vcov aligns with coefs
 vc <- vcov(object)
 
 # Parameter selection
 pnames <- names(coefs)
 if (missing(parm)) {
  parm <- pnames
 } else if (is.numeric(parm)) {
  parm <- pnames[parm]
 }
 
 # Critical value (Standard Normal for Cox/Large sample)
 a <- (1 - level) / 2
 a <- c(a, 1 - a)
 crit <- stats::qnorm(1 - (1 - level) / 2)
 
 # Format column names
 pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
 ci <- array(NA, dim = c(length(parm), 2), dimnames = list(parm, pct))
 
 # Computation
 ses <- sqrt(diag(vc))
 names(ses) <- rownames(vc)
 
 # Safe matching of parameters to SEs
 idx <- match(parm, names(coefs))
 est_sub <- coefs[idx]
 se_sub <- ses[match(parm, names(ses))]
 
 ci[] <- est_sub + se_sub %o% c(-crit, crit)
 
 return(ci)
}

#' Summarizing Cox PH Mixture Fits
#'
#' @description
#' \code{summary} method for class \code{coxphMixture}. Provides a detailed summary 
#' of the fitted model, including coefficients, hazard ratios, standard errors, 
#' z-statistics, and p-values for both the outcome model and the mismatch model.
#'
#' @param object An object of class \code{coxphMixture}.
#' @param conf.int The confidence level for the confidence intervals of the hazard ratios.
#' @param scale Scale factor for the standard errors (default is 1).
#' @param ... Additional arguments.
#'
#' @return An object of class \code{summary.coxphMixture} containing:
#' \item{call}{The function call.}
#' \item{n}{Total number of observations.}
#' \item{nevent}{Number of events.}
#' \item{coefficients}{Matrix of coefficients for the outcome model.}
#' \item{m.coefficients}{Matrix of coefficients for the mismatch model.}
#' \item{conf.int}{Matrix of confidence intervals for the hazard ratios.}
#' \item{logtest}{Log-likelihood information (Outcome Model).}
#' \item{avgcmr}{The average posterior probability of a correct match.}
#'
#' @export
summary.coxphMixture <- function(object, conf.int = 0.95, scale = 1, ...) {
 
 # Outcome Model Coefficients
 beta <- object$coefficients
 n_beta <- length(beta)
 
 # Extract SEs from the top-left block of the variance matrix
 # object$var includes both beta and gamma. 
 # We assume the first n_beta rows/cols correspond to beta.
 se_beta <- sqrt(diag(object$var)[1:n_beta])
 z_beta <- beta / se_beta
 p_beta <- 2 * (1 - stats::pnorm(abs(z_beta)))
 
 coef_mat <- cbind(
  coef = beta,
  `exp(coef)` = exp(beta),
  `se(coef)` = se_beta,
  z = z_beta,
  `Pr(>|z|)` = p_beta
 )
 rownames(coef_mat) <- names(beta)
 
 # Confidence Intervals for Hazard Ratios
 z_crit <- stats::qnorm((1 + conf.int) / 2, 0, 1)
 tmp <- cbind(
  exp(beta),
  exp(-beta),
  exp(beta - z_crit * se_beta),
  exp(beta + z_crit * se_beta)
 )
 colnames(tmp) <- c("exp(coef)", "exp(-coef)", 
                    paste("lower .", round(100 * conf.int, 2), sep = ""),
                    paste("upper .", round(100 * conf.int, 2), sep = ""))
 rownames(tmp) <- names(beta)
 
 # Mismatch Model Coefficients
 gamma <- object$m.coefficients
 m_coef_mat <- NULL
 
 if (length(gamma) > 0) {
  # Extract SEs for gamma (remainder of the diagonal)
  se_gamma <- sqrt(diag(object$var)[(n_beta + 1):nrow(object$var)])
  z_gamma <- gamma / se_gamma
  p_gamma <- 2 * (1 - stats::pnorm(abs(z_gamma)))
  
  m_coef_mat <- cbind(
   Estimate = gamma,
   `Std. Error` = se_gamma,
   `z value` = z_gamma,
   `Pr(>|z|)` = p_gamma
  )
  rownames(m_coef_mat) <- names(gamma)
 }
 
 # Construct Result
 res <- list(
  call = object$call,
  n = object$n,
  nevent = object$nevent,
  coefficients = coef_mat,
  m.coefficients = m_coef_mat,
  conf.int = tmp,
  avgcmr = mean(object$match.prob, na.rm = TRUE),
  iter = length(object$objective)
 )
 
 class(res) <- "summary.coxphMixture"
 return(res)
}

#' @export
print.summary.coxphMixture <- function(x, 
                                       digits = max(3L, getOption("digits") - 3L),
                                       signif.stars = getOption("show.signif.stars"), 
                                       ...) {
 cat("\nCall:\n")
 dput(x$call)
 
 cat("\n--- Outcome Model (Cox PH) ---\n")
 stats::printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, P.values = TRUE, has.Pvalue = TRUE)
 
 cat("\n--- Hazard Ratios & Confidence Intervals ---\n")
 print(x$conf.int, digits = digits)
 
 if (!is.null(x$m.coefficients)) {
  cat("\n--- Mismatch Indicator Model ---\n")
  stats::printCoefmat(x$m.coefficients, digits = digits, signif.stars = signif.stars, P.values = TRUE, has.Pvalue = TRUE)
 }
 
 cat("\nAverage Estimated Correct Match Rate:", format(x$avgcmr, digits = digits), "\n")
 cat("Events:", x$nevent, " / Total:", x$n, "\n")
 cat("Iterations:", x$iter, "\n\n")
 
 invisible(x)
}

#' Print a coxphMixture Object
#'
#' @param x An object of class \code{coxphMixture}.
#' @param digits The number of significant digits to use.
#' @param ... Additional arguments.
#' @export
print.coxphMixture <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
 cat("Call:\n")
 dput(x$call)
 
 cat("\nOutcome Model Coefficients:\n")
 print(format(x$coefficients, digits = digits), quote = FALSE)
 
 if (length(x$m.coefficients) > 0) {
  cat("\nMismatch Model Coefficients:\n")
  print(format(x$m.coefficients, digits = digits), quote = FALSE)
 }
 
 cat("\nLikelihood ratio test (model=outcome) not available due to pseudo-likelihood.\n")
 cat("n=", x$n, ", number of events=", x$nevent, "\n")
 invisible(x)
}

#' Model Predictions for coxphMixture Objects
#'
#' @description
#' Compute fitted values and predictions for the outcome model component of the mixture.
#' The predictions are conditional on the latent status being a "correct match".
#'
#' @param object An object of class \code{coxphMixture}.
#' @param newdata Optional new data frame. If missing, predictions are for the original data.
#' @param type The type of prediction.
#' \itemize{
#'   \item \code{"lp"}: Linear predictor (eta = X * beta).
#'   \item \code{"risk"}: Risk score (exp(eta)).
#'   \item \code{"expected"}: Expected number of events (approximate).
#'   \item \code{"survival"}: Survival probability at the observed times.
#' }
#' @param se.fit Logical; whether to compute standard errors (based on the sandwich/Louis variance).
#' @param na.action Function to handle missing values in \code{newdata}.
#' @param reference Reference for centering (currently ignored, defaults to uncentered).
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' When \code{newdata} is supplied, the function constructs the model matrix using
#' the terms from the original fit. Standard errors are computed using the
#' estimated variance-covariance matrix of the mixture model coefficients.
#' 
#' For \code{type = "expected"} and \code{"survival"}, the function reconstructs the 
#' cumulative baseline hazard step function \eqn{\Lambda_0(t)} using the Breslow 
#' estimator stored in the object and evaluates it at the time points found in 
#' \code{newdata}.
#'
#' @return A vector or matrix of predictions, or a list containing \code{fit} and 
#' \code{se.fit} if standard errors are requested.
#'
#' @importFrom stats terms delete.response model.frame model.matrix predict approxfun
#' @export
predict.coxphMixture <- function(object, newdata,
                                 type = c("lp", "risk", "expected", "survival"),
                                 se.fit = FALSE,
                                 na.action = stats::na.pass,
                                 reference = "strata",
                                 ...) {
 
 type <- match.arg(type)
 
 # Construct Design Matrix (X) and Identify Time (if needed)
 new_time <- NULL
 
 if (missing(newdata)) {
  # Use original X if available
  if (!is.null(object$x)) {
   X <- object$x
  } else if (!is.null(object$model)) {
   # Reconstruct from model frame
   trms <- stats::delete.response(object$terms)
   X <- stats::model.matrix(trms, object$model)
   if (attr(object$terms, "intercept") == 1) X <- X[, -1, drop = FALSE]
  } else {
   stop("Original model matrix not found. Please refit with 'x = TRUE' or provide 'newdata'.")
  }
  
  # For original data, Lambdahat0 is already aligned 1:1 with observations
  base_haz_vals <- object$Lambdahat0
  
 } else {
  # Process New Data
  trms <- stats::delete.response(object$terms)
  mf <- stats::model.frame(trms, data = newdata, na.action = na.action)
  X <- stats::model.matrix(trms, mf)
  if (attr(object$terms, "intercept") == 1) X <- X[, -1, drop = FALSE]
  
  # If type is expected/survival, we need the Time variable from newdata
  if (type %in% c("expected", "survival")) {
   # Without the Surv object, we rely on the formula for the time variable.
   surv_obj_call <- object$call$formula[[2]] # LHS of formula (Surv(time, status))
   
   # Try to extract time variable name. 
   # Usually the 2nd element of Surv() call is time.
   time_var_name <- as.character(surv_obj_call[[2]])
   
   if (time_var_name %in% names(newdata)) {
    new_time <- newdata[[time_var_name]]
   } else {
    stop("Could not identify the time variable in 'newdata' required for expected/survival predictions.")
   }
   
   # Reconstruct Step Function for Baseline Hazard
   # object$Lambdahat0 corresponds to the sorted event times of the original data.
   # We need the training times to build the step function.
   if (is.null(object$y)) {
    stop("The training response 'y' is missing from the object. Please refit with 'y = TRUE' to support expected/survival prediction on new data.")
   }
   
   # Extract training times and hazard
   train_times <- object$y[, "time"]
   
   # We create a step function: H(t). method="constant" gives right-continuous step function (f=0)
   # or appropriate step behavior. Breslow is usually constant between events.
   # We verify if Lambdahat0 is already sorted or 1:1.
   # In 'coxphMixture', Lambdahat0 is expanded to match 'y'. 
   # We reduce it to unique times for the step function.
   
   ord <- order(train_times)
   sorted_times <- train_times[ord]
   sorted_haz <- object$Lambdahat0[ord]
   
   # Create approximation function (Step function)
   # rule=2 means extrapolate constant value outside range (keep max hazard for t > max_t)
   haz_step_fun <- stats::approxfun(sorted_times, sorted_haz, method = "constant", rule = 2)
   
   base_haz_vals <- haz_step_fun(new_time)
  }
 }
 
 # Compute Linear Predictor and Risk
 beta <- object$coefficients
 
 if (ncol(X) != length(beta)) {
  stop(paste("Dimension mismatch: Design matrix has", ncol(X), 
             "columns but model has", length(beta), "coefficients."))
 }
 
 lp <- as.vector(X %*% beta)
 
 if (type == "lp") {
  pred <- lp
 } else if (type == "risk") {
  pred <- exp(lp)
 } else if (type %in% c("expected", "survival")) {
  
  # Cumulative Hazard = Lambda0(t) * exp(lp)
  expected <- base_haz_vals * exp(lp)
  
  if (type == "expected") {
   pred <- expected
  } else {
   pred <- exp(-expected)
  }
 }
 
 # Compute Standard Errors
 if (se.fit) {
  n_beta <- length(beta)
  cov_beta <- object$var[1:n_beta, 1:n_beta, drop = FALSE]
  
  # Variance of LP: diag(X * Sigma * X')
  var_lp <- rowSums((X %*% cov_beta) * X)
  se_lp <- sqrt(var_lp)
  
  if (type == "lp") {
   se <- se_lp
  } else if (type == "risk") {
   # Delta method: Var(exp(lp)) = (exp(lp))^2 * Var(lp)
   se <- pred * se_lp
  } else {
   # Standard errors for expected/survival involve the variance of the baseline hazard,
   # which is computationally expensive in mixture models (requires full Louis matrix integration).
   # We return NA to avoid misleading users with partial SEs.
   warning("Standard errors for 'expected' and 'survival' types are not currently available for mixture models.")
   se <- rep(NA, length(pred))
  }
  
  return(list(fit = pred, se.fit = se))
 }
 
 return(pred)
}