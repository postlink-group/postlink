#' Summarize a "glm_ele" Object
#' @description Summarize results from a `glm_ele()` object
#'
#' @param object the result of a call to `glm_ele()`
#' @param ... for additional summary arguments
#'
#' @returns a list of results from the function called depending on the "family" specified.
#' \item{call}{the matched call}
#' \item{family}{the assumed type of (outcome) regression model}
#' \item{coefficients}{a matrix with the outcome model's coefficient estimates, standard errors, t or z values, and p-values}
#' \item{dispersion}{the dispersion parameter estimate when the family is a Generalized Linear Model}
#'
#' @export
summary.glm_ele <- function(object,...){
 tab <- list()

 for (i in 1:nrow(object$coefficients)){
  if (object$family == "gaussian" | object$family == "gamma"){
  tval <- object$coefficients[i,]/object$standard.errors[i,]
  df.residual <- object$df.residual
  pval <- 2*pt(abs(tval),df=df.residual, lower.tail = FALSE)
  TAB <- cbind(object$coefficients[i,], object$standard.errors[i,],
               tval, pval)
  colnames(TAB) <- c("Estimate","Std. Error", "t value", "Pr(>|t|)")
  rownames(TAB) <- rownames(TAB)
  }
  
  if (object$family == "poisson" | object$family == "binomial"){
   zval <- object$coefficients[i,]/object$standard.errors[i,]
   pval <- 2 * (1 - pnorm(abs(zval)))
   
   TAB <- cbind(object$coefficients[i,], object$standard.errors[i,],
                zval, pval)
   colnames(TAB) <- c("Estimate","Std. Error", "z value", "Pr(>|z|)")
   rownames(TAB) <- rownames(TAB)
  }
  tab[[i]] <- TAB
 }
 
 names(tab) <- rownames(object$coefficients)
 
 vals <- list(call = object$call, family = object$family,
                coefficients = tab)
 
 if (object$family %in% c("gamma", "gaussian")){
  vals$dispersion <- object$dispersion
 }

 class(vals)    <- "summary.glm_ele"
 vals
}