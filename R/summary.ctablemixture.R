#' Summarize a `ctable_mixture` Object
#' @description
#' Summarize the results in the `ctable_mixture` object.
#'
#' @param object an object of class "ctablemixture".
#' @param ... for additional `summary` arguments.
#'
#' @returns a list of results.
#' \item{call}{the matched call.}
#' \item{phat}{a matrix of the estimated cell probabilities.}
#' \item{vcov_phat}{the variance-covariance matrix of phat.}
#' \item{ftable}{a flat contingency table in matrix format containing the estimated counts
#' of each combination of the variables’ levels.}
#' \item{chi2_statistic}{the Pearson’s chi-squared test statistic based on
#' the `ftable` value}
#'
#' @export
summary.ctablemixture <- function(object,...){

 cat("Call:\n")
 print(object$call, quote = F)
 cat("\n")

 cat("Estimated Cell Counts:\n")
 print(object$ftable)
 cat("\n")

 chi2 <- chisq.test(fit$ftable)

 cat("Pearson's Chi-Squared Test:\n")
 cat("X-squared = ", chi2$statistic, ",", sep = "")
 cat(" df = ", chi2$parameter, ",", sep = "")

 chi2pval <- chi2$p.value
 if(chi2pval < 2.2e-16){
 cat(" p-value < 2.2e-16")
 } else{
  cat(" p-value =", chi2$statistic)
 }
 cat("\n")

 object <- list(call = object$call, phat = object$phat,
                vcov_phat = object$vcov_phat, ftable = object$ftable,
                chi2_statistic = chi2$statistic)

 invisible(object)
}
