#' Print method for summary.glm_mixtureBayes
#'
#' @param x An object of class \code{"summary.glm_mixtureBayes"}.
#' @param digits Significant digits to use in printing.
#' @param ... Additional arguments (unused).
#' @return The input \code{x}, invisibly.
#' @export
print.summary.glm_mixtureBayes <- function(x, digits = max(3L, getOption("digits") - 3L),
                                      signif.stars = getOption("show.signif.stars"),...){
  cat("Call:", sep="\n")
  print(x$call,quote=F)
  cat(" ", sep="\n")
  cat("Family:", x$family, " ", sep="\n")
  
  cat("(For Correct Matches):", sep="\n")
  
  cat("Outcome Model Coefficients:", sep="\n")
  printCoefmat(x$coefficients,quote=F, digits = digits,
               signif.stars = signif.stars)
  cat(" ", sep="\n")
  
  if (x$family %in% c("gamma", "gaussian")){
    cat("Dispersion:", sep="\n")
    print(format(signif(x$dispersion, digits)), print.gap = 2, quote = F)
    cat("\n")
  }
  
  cat("(For Mismatches):", sep="\n")
  
  cat("Outcome Model Coefficients:", sep="\n")
  printCoefmat(x$m.coefficients,quote=F, digits = digits,
               signif.stars = signif.stars)
  cat(" ", sep="\n")
  
  if (x$family %in% c("gamma", "gaussian")){
    cat("Dispersion:", sep="\n")
    print(format(signif(x$m.dispersion, digits)), print.gap = 2, quote = F)
    cat("\n")
  }
  
  invisible(x)
}