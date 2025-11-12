#' Print a `ctable_mixture` Object
#' @description
#' Print the call and input data for the `ctable_mixture` object.
#'
#' @param object an object of class "ctablemixture".
#' @param digits the number of significant digits to use for printing.
#' @param ... for additional `print` arguments.
#'
#' @returns Invisibly returns the object from `ctable_mixture()`
#' that is provided as an argument.
print.ctablemixture <- function(object, digits = getOption("digits"),...){
 cat("Call: ")
 print(object$call, quote = F, digits = digits)

 cat("Mismatch Rate:", object$m.rate)
 cat("\n")
 cat("\n")

 cat("Number of Samples:", sum(object$ftable))
 cat("\n")

 cat("Table Dimensions:", dim(object$ftable)[1], "rows", "x", dim(object$ftable)[2], "columns")
 cat("\n")

 invisible(object)
}
