#' Print Method for `adjELE` Objects
#'
#' Provides a concise summary of the adjustment object created by \code{\link{adjELE}},
#' including linkage error assumptions, blocking structure, and weight estimation settings.
#'
#' @param x An object of class \code{adjELE}.
#' @param digits Integer; the number of significant digits to use when printing
#'   numeric values. Defaults to 3.
#' @param ... Additional arguments passed to methods.
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @details
#' This method inspects the internal structure of the adjustment object.
#' It calculates summaries for mismatch rates and audit sizes (e.g., means/ranges)
#' if they vary across blocks, providing a snapshot of the error assumption complexity.
#' It safely handles cases where the reference data is missing or empty.
#'
#' @examples
#' # Example: Using the included brfss demonstration dataset
#' data(brfss, package = "postlink")
#'
#' adj_object <- adjELE(linked.data = brfss,
#'                     m.rate = unique(brfss$m.rate),
#'                     blocks = imonth,
#'                     weight.matrix = "BLUE")
#' print(adj_object)
#'
#' @export
print.adjELE <- function(x, digits = 3, ...) {

 # 1. Header (ASCII for CRAN compatibility)
 cat("\n--- Adjustment Object: Exchangeable Linkage Errors (Chambers, 2009) ---\n")

 # 2. Data Summary
 has_data <- FALSE
 n_obs <- 0

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
  cat("\n    Observations:  ", format(n_obs, big.mark = ","))
  cat("\n    Storage:        Reference (Environment)\n")
 } else {
  cat("\n    Status:         Not available / Empty\n")
 }

 # 3. Method Specification
 cat("\n* Specification:")

 # Weight Matrix Method
 w_mat <- if (!is.null(x$weight.matrix)) x$weight.matrix else "Unknown"
 cat("\n    Weight Matrix:  ", w_mat)

 # Blocking Structure
 # Calculate unique blocks safely (ignoring NAs)
 if (!is.null(x$blocks)) {
  n_unique_blocks <- length(unique(x$blocks[!is.na(x$blocks)]))
  if (n_unique_blocks == 1) {
   cat("\n    Blocks:         Single Block (Global assumption)")
  } else {
   cat("\n    Blocks:         ", format(n_unique_blocks, big.mark = ","), " distinct blocks", sep = "")
  }
 } else {
  cat("\n    Blocks:         None specified")
 }

 # 4. Mismatch Rates Summary
 cat("\n    Mismatch Rate:  ")
 if (!is.null(x$m.rate)) {
  rates <- x$m.rate
  # Check if constant
  if (length(unique(rates)) == 1) {
   cat(format(rates[1], digits = digits), "(Constant)")
  } else {
   # Summary for variable rates
   mean_r <- mean(rates, na.rm = TRUE)
   min_r <- min(rates, na.rm = TRUE)
   max_r <- max(rates, na.rm = TRUE)
   cat(sprintf("Variable (Mean: %.*f, Range: %.*f - %.*f)",
               digits, mean_r, digits, min_r, digits, max_r))
  }
 } else {
  cat("None specified")
 }

 # 5. Audit Size Summary
 if (!is.null(x$audit.size)) {
  audits <- x$audit.size
  total_audit <- sum(unique(audits), na.rm = TRUE)

  # If audit.size is length 1 or constant
  if (length(unique(audits)) == 1) {
   # we report the Total if computable, or the raw value if global.
   cat("\n    Audit Sample:   Global size", format(audits[1], big.mark = ","))
  } else {
   # If variable, likely specific per block.
   # It's safer to describe variability than try to guess 'Total' without block IDs here
   cat("\n    Audit Sample:   Variable (Range:", min(audits), "-", max(audits), ")")
  }
 } else {
  cat("\n    Audit Sample:   None (Using known rates)")
 }

 cat("\n\n")
 invisible(x)
}
