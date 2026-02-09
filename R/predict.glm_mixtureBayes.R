#' Predictions from a glm_mixtureBayes model
#'
#' @param object A \code{glm_mixtureBayes} model object.
#' @param newdata A data frame of new observations for which to predict.
#' @param type Either \code{"link"} or \code{"response"}, indicating the scale of predictions (linear predictor or original response scale).
#' @param se.fit Logical; if \code{TRUE}, also return standard errors of predictions.
#' @param interval Either \code{"none"} or \code{"credible"}, indicating whether to compute a credible interval for predictions.
#' @param level Probability level for the credible interval (default 0.95).
#' @param na.action How to handle missing values in \code{newdata}. Default is \code{na.pass}.
#' @param ... Not used.
#' @return If \code{se.fit = FALSE} and \code{interval = "none"}, a numeric vector of predicted values. If \code{se.fit = TRUE} or \code{interval = "credible"}, a matrix of predictions with columns for the fit, standard error, and lower/upper bounds of the credible interval.
#' @export
predict.glm_mixtureBayes <- function(object, newdata,
                                type = c("link", "response"),
                                se.fit = FALSE,
                                interval = c("none", "credible"), level = 0.95,
                                na.action = stats::na.pass, ...){
  type <- match.arg(type)
  interval <- match.arg(interval)

  # have not incorporated scenario when newdata is NULL
  X_data <- stats::model.matrix(stats::as.formula(object$call$formula), data = newdata)

  mean_coef <- apply(object$estimates$coefficients, 2, mean)

  # obtain predictions
  # have not incorporated na.action
  if(type == "link"){
    predictions <- X_data %*% mean_coef
  }

  if(type == "response"){
    predictions <- switch(object$family,
                          gaussian = X_data %*% mean_coef, # identity link
                          gamma = exp(X_data %*% mean_coef), # log link
                          poisson = exp(X_data %*% mean_coef), # log link
                          binomial = stats::plogis(X_data %*% mean_coef)) # logit link
  }

  if(se.fit == FALSE & interval == "none"){
    return(c(predictions))
  }

  # obtain se.fit
  if(type == "link"){
    # rows are new data points, columns are draws
    all_predictions <- X_data %*% t(object$estimates$coefficients)
  }

  if(type == "response"){
    # rows are new data points, columns are draws
    all_predictions <- switch(object$family,
                       gaussian = X_data %*% t(object$estimates$coefficients),
                       gamma = exp(X_data %*% t(object$estimates$coefficients)),
                       poisson = exp(X_data %*% t(object$estimates$coefficients)), # log link
                       binomial = stats::plogis(X_data %*% t(object$estimates$coefficients))) # logit link
  }

  se.predictions <- apply(all_predictions, 1, stats::sd)

  # obtain confidence interval
  if(interval == "none"){
    vals <- cbind(fit = predictions, se.fit = se.predictions)
    colnames(vals) <- c("fit", "se.fit")
    return(vals)
  }

  alpha <- 1 - level
  ci <- t(apply(all_predictions, 1, stats::quantile, probs = c(alpha/2, 1-alpha/2)))

  lower <- ci[,1]
  upper <- ci[,2]

  if(se.fit == FALSE & interval == "credible"){
    vals <- cbind(fit = predictions, lower = lower, upper = upper)
    colnames(vals) <- c("fit", paste(alpha/2*100, "%"),
                             paste((1 - alpha/2)*100, "%"))
    return(vals)
  }

  if(se.fit == TRUE & interval == "credible"){
    vals <- cbind(fit = predictions, se.fit = se.predictions,
                 lower = lower, upper = upper)
    colnames(vals) <- c("fit", "se.fit", paste(alpha/2*100, "%"),
                          paste((1 - alpha/2)*100, "%"))
    return(vals)
  }
}
