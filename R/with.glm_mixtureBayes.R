#' Multiple Imputation Pooling from Linkage Posterior Samples
#'
#' @description
#' Pool regression estimates across many "completed datasets" that arise when
#' each posterior draw of record-linkage indicators \eqn{z} selects a subset of
#' records to treat as truly linked. For each retained draw, the function fits
#' the requested model (LM/GLM) on the linked subset, collects coefficients and
#' their covariance, and combines them using Rubin's rules (vector/matrix form).
#'
#' @param object A `glm_mixtureBayes` model object containing posterior linkage samples.
#' @param data A data.frame with all candidate records in the same row order as used in the model.
#' @param formula Model formula (defaults to `object$formula`).
#' @param family GLM family (defaults to `object$family`).
#' @param min_n Minimum sample size required to fit the model for a given draw.
#'   Defaults to \code{p + 1}, where \code{p} is the number of columns in the model matrix.
#' @param quietly If \code{TRUE}, suppress errors from individual failed fits and skip them.
#'
#' @return An object of class \code{"mi_link_pool"} with components:
#' \itemize{
#'   \item \code{m}: number of retained imputations (successful fits).
#'   \item \code{coef}: pooled coefficient estimates (vector, length \eqn{p}).
#'   \item \code{vcov}: total covariance matrix \eqn{T = \bar U + (1 + 1/m) B} (\eqn{p \times p}).
#'   \item \code{se}: pooled standard errors (diagonal of \code{vcov}).
#'   \item \code{ci95}: 95\% confidence intervals (matrix with columns \code{lwr}, \code{upr}).
#'   \item \code{Ubar}: average within-imputation covariance (\eqn{\bar U}).
#'   \item \code{B}: between-imputation covariance of point estimates.
#'   \item \code{lambda}: fraction of missing information per parameter.
#'   \item \code{df}: per-parameter degrees of freedom (Rubin 1987 approximation).
#'   \item \code{kept_draws}: indices of posterior draws that were used.
#' }
#'
#' @details
#' Draws yielding no linked records or too few observations are **skipped** rather
#' than filled with zeros—this avoids biasing pooled estimates. The function pools
#' all parameters jointly (vectorized Rubin's rules), preserving cross-parameter
#' covariance, which is important for contrasts and multivariate inference.
#'
#' @examples
#' \dontrun{
#' # Suppose z_mat is G x n, data is a data.frame, and you fit y ~ x with Gaussian errors:
#' res <- with.glm_mixtureBayes(object,
#'                              data   = dat_linked,
#'                              formula = y ~ x,
#'                              family  = gaussian())
#' print(res)
#' coef(res)
#' sqrt(diag(res$vcov))
#' }
#'
#' @export
#' @importFrom stats coef vcov cov lm glm qt model.frame model.matrix

with.glm_mixtureBayes <- function(object, data, formula = object,
                                  family = object$family, min_n = NULL, quietly = TRUE) {
  # Extract the matrix of linkage indicators from the fitted model:
  z_samples <- object$m_samples

  # ---- Helper: extract a length-n 0/1 vector for the g-th draw ----------------
  collapse_z_g <- function(z_samples, g){
    if (is.list(z_samples)){
      # Each element should be a matrix-like object with G rows.
      parts <- lapply(z_samples, function(Z){
        if (is.null(dim(Z)))
          stop("Each list element of `z_samples` must be a matrix-like object.")
        Z[g, ]
      })
      as.integer(do.call(c, parts))
    } else {
      # Matrix case: G x n
      if (is.null(dim(z_samples)) || length(dim(z_samples)) != 2L)
        stop("`z_samples` must be a G x n matrix, or a list of such matrices.")
      as.integer(z_samples[g, ])
    }
  }

  if (missing(formula) || is.null(formula)) {
    if (!is.null(object$call$formula)) {
      formula <- eval(object$call$formula, envir = parent.frame())
    } else {
      stop("Formula not found: please provide it explicitly or ensure object$call$formula exists.")
    }
  }

  # ---- Helper: fit once (LM for gaussian, GLM otherwise); return coef & vcov ---
  fit_once <- function(df) {
    if (inherits(family, "family") && family$family == "gaussian"){
      fit <- lm(formula, data = df)
    } else {
      fit <- glm(formula, data = df, family = family)
    }
    list(coef = stats::coef(fit), vcov = stats::vcov(fit))
  }

  # ---- Determine p and default min_n ------------------------------------------
  mf <- model.frame(formula, data = data)
  Xtmp <- model.matrix(formula, mf)
  p <- ncol(Xtmp)
  if (is.null(min_n)) min_n <- p + 1L

  # ---- Number of posterior draws G --------------------------------------------
  G <- if (is.list(z_samples)) nrow(z_samples[[1]]) else nrow(z_samples)
  if (length(G) == 0L || is.null(G) || !is.finite(G))
    stop("Unable to infer the number of posterior draws (G).")

  coefs_list <- list()
  vcovs_list <- list()
  kept <- logical(G)

  # ---- Iterate over posterior draws; fit on linked subset if size >= min_n ----
  for (g in seq_len(G)){
    zg <- collapse_z_g(z_samples, g)
    if (!all(zg %in% c(1L, 2L)))
      stop("`z_samples` must contain only 1/2 indicators.")

    if (sum(zg) < min_n) next
    dat_g <- data[which(zg == 1L), , drop = FALSE]

    # Attempt a fit; skip if it fails (e.g., perfect separation in GLM, collinearity)
    res <- try(fit_once(dat_g), silent = quietly)
    if (!inherits(res, "try-error")){
      coefs_list[[length(coefs_list)+1L]] <- res$coef
      vcovs_list[[length(vcovs_list)+1L]] <- res$vcov
      kept[g] <- TRUE
    }
  }

  m <- length(coefs_list)
  if (m == 0L)
    stop("No valid imputations: all draws failed or had too few linked records.")

  # ---- Align parameter names; build m x p coefficient matrix -------------------
  parnames <- names(coefs_list[[1L]])
  coefs_mat <- do.call(rbind, lapply(coefs_list, function(b) b[parnames]))

  # ---- Within-imputation covariance (average), and between-imputation cov -----
  Ubar <- Reduce("+", vcovs_list) / m           # average within-imputation covariance
  B    <- stats::cov(coefs_mat)                 # between-imputation covariance of estimates

  # ---- Rubin's rules (vectorized): total covariance T = Ubar + (1 + 1/m) * B ---
  Tmat <- Ubar + (1 + 1/m) * B

  # ---- Pooled estimates, standard errors, fraction missing information --------
  qbar   <- colMeans(coefs_mat)                 # pooled point estimates
  se     <- sqrt(diag(Tmat))                    # pooled standard errors
  lambda <- diag((1 + 1/m) * B) / diag(Tmat)    # fraction of missing information

  # ---- Degrees of freedom (Rubin 1987 "old" df, per-parameter) -----------------
  # r = ((1 + 1/m) * B) / Ubar, computed on diagonals (per parameter)
  r     <- diag((1 + 1/m) * B) / diag(Ubar)
  dfold <- (m - 1) * (1 + 1/r)^2
  # Conservative floor for numerical stability in small m:
  df    <- pmax(dfold, 3)

  # ---- 95% t-intervals ---------------------------------------------------------
  tcrit <- stats::qt(0.975, df = df)
  lwr   <- qbar - tcrit * se
  upr   <- qbar + tcrit * se
  ci95  <- cbind(lwr = lwr, upr = upr)
  rownames(ci95) <- names(qbar)

  out <- list(
    m          = m,
    coef       = qbar,
    vcov       = Tmat,
    se         = se,
    ci95       = ci95,
    Ubar       = Ubar,
    B          = B,
    lambda     = lambda,
    df         = df,
    kept_draws = which(kept)
  )
  class(out) <- c("mi_link_pool_glm", "mi_link_pool")
  return(out)
}

#' @export
print.mi_link_pool_glm <- function(x, digits = max(3L, getOption("digits") - 2L), ...){
  # Compact, human-friendly summary for console output.
  cat("Multiple Imputation (from linkage draws):\n")
  cat("  Retained imputations (m):", x$m, "\n\n")

  est <- x$coef
  se  <- x$se
  ci  <- x$ci95
  df  <- x$df
  tab <- cbind(Estimate = est, Std.Error = se,
               `CI.lwr` = ci[, "lwr"], `CI.upr` = ci[, "upr"],
               df = df)
  print(round(tab, digits))
  invisible(x)
}
