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
#' where the linked data is unspecified (NULL). It also prints the user-specified
#' priors or outlines the defaults that will be used.
#'
#' @examples
#' data(lifem)
#'
#' # lifem data preprocessing
#' # For computational efficiency in the example, we work with a subset of the lifem data.
#' lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
#' lifem_small <- rbind(
#'   head(subset(lifem, hndlnk == 1), 100),
#'   head(subset(lifem, hndlnk == 0), 20)
#' )
#'
#' adj_obj <- adjMixBayes(
#'   linked.data = lifem_small,
#'   priors = list(theta = "beta(2, 2)")
#' )
#'
#' # Implicitly calls print.adjMixBayes()
#' adj_obj
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
 } else {
  cat("\n    Status:       None specified (NULL)\n")
 }

 cat("\n* Priors:")
 if (!is.null(x$priors) && length(x$priors) > 0) {
  cat("\n    User-specified overrides:\n")
  for (p in names(x$priors)) {
   cat(sprintf("      %-10s : %s\n", p, x$priors[[p]]))
  }
  cat("    (Unspecified parameters will use defaults below)\n")
 } else {
  cat("\n    Status:       None specified. Using symmetric defaults.\n")
 }

 cat("    Defaults applied during fitting:\n")
 cat("      GLM:        beta ~ normal(0,5) [binomial: normal(0,2.5)]\n")
 cat("      Survival:   beta ~ normal(0,5) [weibull: normal(0,2)]\n")
 cat("      Mix Weight: theta ~ beta(1,1)\n")

 cat("\n")
 invisible(x)
}
