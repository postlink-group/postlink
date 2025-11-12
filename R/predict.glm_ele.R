#' Predictions From a `glm_ele` Object
#' @description Obtain predictions from a `glm_ele()` object.
#'
#' @param object the result of a call to `glm_ele()`
#' @param newdata optional new data to obtain predictions for. The original data is used by default.
#' @param type the type of prediction. The choices are predictions on the scale
#' of the linear predictors ("link") or response ("response").
#' @param ... for future predict arguments
#'
#' @returns a vector of predictions based on arguments specified.
#'
#' @export
predict.glm_ele <- function(object, newdata, type = c("response", "link")){
 if (!missing(newdata)){
  coef <- object$coefficients
  X <- model.matrix(~ ., data = newdata)
  lp <- X %*% coef
 } else{
  lp <- object$linear.predictors
 }

 colnames(lp) <- rownames(object$coefficients)
 type <- match.arg(type)
 if(type == "link"){
  return(lp)
 }

 if(type == "response"){
  vals <- switch(object$family,
         gaussian = lp,
         poisson = exp(lp),
         gamma = exp(lp),
         binomial = exp(lp) / (1 + exp(lp)))
  return(vals)
 }
}
