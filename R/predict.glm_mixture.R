#' Predictions from a `glm_mixture` Object
#' @description Obtain predictions from a `glm_mixture` object.
#'
#' @param object the result of a call to `glm_mixture()`
#' @param newdata optional new data to obtain predictions for. The original data is used by default.
#' @param type the type of prediction. The choices are predictions on the scale
#' of the linear predictors ("link") or response ("response").
#' @param se.fit an indicator of whether standard errors are required (TRUE or FALSE)
#' @param dispersion the dispersion to be assumed when computing standard errors.
#' By default, it is the value returned in the fit object.
#' @param interval type of interval calculation.
#' @param level required confidence level if interval is set to "confidence".
#' @param na.action a function for what to do with missing values in `newdata`.
#' The default is to predict "NA".
#' @param ... for future predict arguments
#'
#' @returns a vector or matrix of predictions based on arguments specified.
#'
#' @examples
#' ## commonness score of first and last names used for linkage
#' mformula <- ~commf + comml
#' ## hand-linked records are considered "safe" matches
#' safematches <- ifelse(lifem$hndlnk =="Hand-Linked At Some Level", TRUE, FALSE)
#' ## overall mismatch rate in the data set is assumed to be ~ 0.05
#' mrate <- 0.05
#' fit <- glm_mixture(age_at_death ~ poly(unit_yob, 3, raw = TRUE), data = lifem,
#'                    family = "gaussian", mformula, safematches, mrate)
#'
#' predict(fit)
#'
#' @export
predict.glm_mixture <- function(object, newdata = NULL,
                                type = c("link", "response"),
                                se.fit = FALSE, dispersion = NULL,
                                interval = c("none", "confidence"), level = 0.95,
                                na.action = na.pass, ...){
  type <- match.arg(type)
  interval <- match.arg(interval)

  # obtain predictions
  glm_object <- object$wfit
  predictions <- unname(predict(object = glm_object, newdata = newdata,
                                type = type, na.action = na.action))

  if(se.fit == FALSE & interval == "none"){
    return(predictions)
  }

  # obtain se.fit
  X_data <- if(is.null(newdata)){model.matrix(glm_object)} else{
    model.matrix(glm_object$formula, data = newdata)
  }

  ncoef <- length(object$coefficients)
  if(is.null(dispersion)){dispersion <- object$dispersion}
  cov_matrix_unscaled <- object$vcov[1:ncoef, 1:ncoef] / dispersion

  var_unscaled <- numeric(nrow(X_data))
  for (i in 1:nrow(X_data)) {
    var_unscaled[i] <- X_data[i, , drop = FALSE] %*% cov_matrix_unscaled %*%
      t(X_data[i, , drop = FALSE])
  }
  se.predictions <- sqrt(dispersion*var_unscaled)

  # obtain confidence interval
  if(interval == "none"){
    return(list(fit = predictions, se.fit = se.predictions,
                residual.scale = sqrt(dispersion),
                df.residual = glm_object$df.residual))
  }

  alpha <- 1 - level
  if(object$family %in% c("gaussian", "gamma")){
    cval <- qt(1-alpha/2, glm_object$df.residual)
  } else{
    cval <- qnorm(1 - alpha/2)
  }

  lower <- predictions - cval * se.predictions
  upper <- predictions + cval * se.predictions

  if(type == "response"){
    lower <- glm_object$family$linkinv(lower)
    upper <- glm_object$family$linkinv(upper)
  }

  if(se.fit == FALSE & interval == "confidence"){
    vals <- cbind(fit = predictions, lower = lower, upper = upper)
    colnames(vals)[2:3] <- c(paste(alpha/2*100, "%"),
                          paste((1 - alpha/2)*100, "%"))
    return(vals)
  }

  if(se.fit == TRUE & interval == "confidence"){
    vals <- list(fit = predictions, se.fit = se.predictions,
                 residual.scale = sqrt(dispersion),
                 df.residual = glm_object$df.residual,
                 lower = lower, upper = upper)
    names(vals)[5:6] <- c(paste(alpha/2*100, "%"),
                          paste((1 - alpha/2)*100, "%"))
    return(vals)
  }

}
