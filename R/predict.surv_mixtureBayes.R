#' Predictions from a surv_mixtureBayes model
#'
#' @param object A \code{surv_mixtureBayes} model object.
#' @param newdata A data frame of new observations for which to predict.
#' @param type Either \code{"link"} or \code{"response"}, indicating the scale of the predictions (linear predictor or original outcome scale).
#' @param se.fit Logical; if \code{TRUE}, also return standard errors of predictions.
#' @param interval Either \code{"none"} or \code{"credible"}, indicating whether to compute a credible interval for predictions.
#' @param level Probability level for the credible interval (default 0.95).
#' @param na.action ... Default is \code{stats::na.pass}.
#' @param ... Not used.
#' @return If \code{se.fit = FALSE} and \code{interval = "none"}, returns a numeric vector of predicted values. If \code{se.fit = TRUE} or \code{interval = "credible"}, returns a matrix with columns for the fit, standard error, and lower/upper bounds of the credible interval.
#' @export
predict.surv_mixtureBayes <- function(object, newdata,
                                      type = c("link", "response"),
                                      se.fit = FALSE,
                                      interval = c("none", "credible"), level = 0.95,
                                      na.action = stats::na.pass, ...) {
 type <- match.arg(type)
 interval <- match.arg(interval)

 # Construct design matrix from formula and new data
 X_data <- stats::model.matrix(stats::as.formula(object$call$formula), data = newdata)

 # Posterior mean of component 1 coefficients
 mean_coef <- apply(object$estimates$coefficients, 2, mean)

 # Compute predictions on link or response scale
 if (type == "link") {
  predictions <- X_data %*% mean_coef
 }

 if (type == "response") {
  predictions <- switch(object$family,
                        gamma = exp(X_data %*% mean_coef),  # log-link mean
                        weibull = {
                         shape <- mean(object$estimates$shape)
                         exp(X_data %*% mean_coef) * gamma(1 + 1 / shape)
                        },
                        stop("Unsupported family for survival mixture model."))
 }

 if (!se.fit && interval == "none") {
  return(c(predictions))
 }

 # Matrix of predictions across posterior draws
 if (type == "link") {
  all_preds <- X_data %*% t(object$estimates$coefficients)
 }

 if (type == "response") {
  all_preds <- switch(object$family,
                      gamma = exp(X_data %*% t(object$estimates$coefficients)),
                      weibull = {
                       shape_draws <- object$estimates$shape
                       eta <- X_data %*% t(object$estimates$coefficients)  # n x S
                       t(apply(eta, 1, function(row, shape) {
                        exp(row) * gamma(1 + 1 / shape)
                       }, shape = shape_draws))
                      },
                      stop("Unsupported family for survival mixture model."))
 }

 # Standard errors
 se.pred <- apply(all_preds, 1, stats::sd)

 if (interval == "none") {
  vals <- cbind(fit = predictions, se.fit = se.pred)
  colnames(vals) <- c("fit", "se.fit")
  return(vals)
 }

 # Credible intervals
 alpha <- 1 - level
 ci <- t(apply(all_preds, 1, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2)))
 lower <- ci[, 1]
 upper <- ci[, 2]

 if (!se.fit && interval == "credible") {
  vals <- cbind(fit = predictions, lower = lower, upper = upper)
  colnames(vals) <- c("fit", paste0(alpha/2 * 100, "%"), paste0((1 - alpha/2) * 100, "%"))
  return(vals)
 }

 if (se.fit && interval == "credible") {
  vals <- cbind(fit = predictions, se.fit = se.pred, lower = lower, upper = upper)
  colnames(vals) <- c("fit", "se.fit", paste0(alpha/2 * 100, "%"), paste0((1 - alpha/2) * 100, "%"))
  return(vals)
 }
}
