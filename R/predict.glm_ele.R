#' Predictions From a `glm_ele` Object
#' @description Obtain predictions from a `glm_ele()` object.
#'
#' @param object the result of a call to `glm_ele()`
#' @param newdata optional new data to obtain predictions for. The original data is used by default.
#' @param weight.matrix the weight matrix of interest if multiple were used ("ratio-type",
#' "Lahiri-Larsen", "BLUE"). By default, it is the one corresponding to the first row of the coefficients matrix.
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
#' @export
predict.glm_ele <- function(object, newdata = NULL, weight.matrix,
                            type = c("link", "response"),
                            se.fit = FALSE, dispersion = NULL,
                            interval = c("none", "confidence"), level = 0.95,
                            na.action = na.pass, ...){
 type <- match.arg(type)
 interval <- match.arg(interval)

 if(missing(weight.matrix)){
  weight.matrix <- rownames(object$coefficients)[1]
 }
 if(length(weight.matrix) != 1){
  stop(("Error: weight.matrix must be of length 1"))
 }
 if(!(weight.matrix %in% rownames(object$coefficients))){
  stop(("Error: weight.matrix was not used for the fitted object"))
 }

 coef <- t(object$coefficients)[, weight.matrix, drop = F]
 vcov_object <- object$vcov[[weight.matrix]]

 if(is.null(newdata)){
  formula <- as.formula(object$call$formula)
  X_data <- model.matrix(formula, data = object$model)
 } else{
  formula <- as.formula(object$call$formula)
  mf <- model.frame(formula, data = newdata, na.action = na.action)
  X_data <- model.matrix(formula, data = mf)
 }

 lp <- X_data %*% coef
 colnames(lp) <- weight.matrix
 # also can be obtained via object$linear.predictors

 # obtain predictions
 if(type == "link"){
  predictions <- lp
 }
 if(type == "response"){
  predictions <- switch(object$family,
                        gaussian = lp,
                        poisson = exp(lp),
                        gamma = exp(lp),
                        binomial = exp(lp) / (1 + exp(lp)))
 }

 if(se.fit == FALSE & interval == "none"){
  return(predictions)
 }

 # obtain se.fit
 ncoef <- ncol(object$coefficients)
 if(is.null(dispersion)){dispersion <- object$dispersion[weight.matrix]}
 cov_matrix_unscaled <- vcov_object[1:ncoef, 1:ncoef] / dispersion

 var_unscaled <- numeric(nrow(X_data))
 for (i in 1:nrow(X_data)) {
  var_unscaled[i] <- X_data[i, , drop = FALSE] %*% cov_matrix_unscaled %*%
   t(X_data[i, , drop = FALSE])
 }
 se.predictions <- sqrt(dispersion*var_unscaled)

 if(interval == "none"){
  return(list(fit = predictions, se.fit = se.predictions,
              residual.scale = sqrt(dispersion),
              df.residual = object$df.residual))
 }

 # obtain confidence interval
 alpha <- 1 - level
 if(object$family %in% c("gaussian", "gamma")){
  cval <- qt(1-alpha/2, object$df.residual)
 } else{
  cval <- qnorm(1 - alpha/2)
 }

 lower <- predictions - cval * se.predictions
 upper <- predictions + cval * se.predictions

 if(type == "response"){
  lower <- switch(object$family,
                  gaussian = lower,
                  poisson = exp(lower),
                  gamma = exp(lower),
                  binomial = exp(lower) / (1 + exp(lower)))
  upper <- switch(object$family,
                  gaussian = upper,
                  poisson = exp(upper),
                  gamma = exp(upper),
                  binomial = exp(upper) / (1 + exp(upper)))
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
               df.residual = object$df.residual,
               lower = lower, upper = upper)
  names(vals)[5:6] <- c(paste(alpha/2*100, "%"),
                        paste((1 - alpha/2)*100, "%"))
  return(vals)
 }
}
