#' Print Method for `adjMixBayes` Objects
#'
#' Provides a concise summary of the Bayesian adjustment object created by
#' \code{\link{adjMixBayes}}.
#'
#' @param x An object of class \code{adjMixBayes}.
#' @param ... Additional arguments passed to methods.
#' @return Invisibly returns \code{x}.
#'
#' @details
#' This method inspects the reference-based data environment to report the number
#' of linked records without copying the full dataset. It safely handles cases
#' where the linked data is unspecified (NULL).
#'
#' @examples
#' # Setup example data
#' data(brfss, package = "postlink")
#'
#' adj_object <- adjMixBayes(linked.data = brfss)
#' print(adj_object)
#'
#' @export
print.adjMixBayes <- function(x, ...) {

 # Header: ASCII used for maximum CRAN compatibility
 cat("\n--- Adjustment Object: Bayesian Mixture (Gutman et al., 2016) ---\n")

 # Check data status safely
 has_data <- FALSE
 n_obs <- 0

 # Verify environment exists and contains non-NULL data
 if (!is.null(x$data_ref) && is.environment(x$data_ref)) {
  if (exists("data", envir = x$data_ref, inherits = FALSE)) {
   stored_data <- x$data_ref$data
   if (!is.null(stored_data) && is.data.frame(stored_data)) {
    has_data <- TRUE
    n_obs <- nrow(stored_data)
   }
  }
 }

 cat("\n* Linked Data:")
 if (has_data) {
  cat("\n    Observations:", format(n_obs, big.mark = ","))
  cat("\n    Storage:      Reference (Environment)\n")
 } else {
  cat("\n    Status:       None specified (NULL)\n")
 }

 cat("\n")
 invisible(x)
}
