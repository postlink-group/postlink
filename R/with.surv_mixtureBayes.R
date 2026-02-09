#' Multiple-imputation pooling for survival mixture models
#'
#' This function implements Rubin's rules to pool regression estimates across
#' multiple draws (imputations) from a Bayesian mixture model with survival outcome.
#' For each draw of the latent match labels (1 = true match), a parametric survival
#' regression (\code{\link[survival]{survreg}}) is fit on the subset of data where
#' the label is 1. Estimates from successful fits are combined across draws.
#'
#' @note \strong{Gamma family is not supported.} If \code{dist = "gamma"}, the function will stop with an error.
#' This is because \code{\link[survival]{survreg}} does not support the gamma distribution.
#' You may instead approximate with \code{"weibull"} or fit your model using Stan.
#'
#' @param object A fitted mixtureBayes model object containing MCMC samples of labels.
#' @param data A data frame containing the variables in the survival model.
#' @param formula A \code{\link{formula}} for the survival model (default: \code{object$call$formula}).
#' @param dist Distribution argument passed to \code{\link[survival]{survreg}} (default from \code{object$family}).
#' @param min_n Minimum number of observations (in subset) required to fit the model (default: \eqn{p+1}, where \eqn{p} is the number of parameters).
#' @param ... Additional arguments passed to \code{\link[survival]{survreg}}.
#'
#' @return An object of class \code{"mi_link_pool"}, a list with components:
#'   \item{call}{The mixture model call.}
#'   \item{pooled}{A data frame of pooled estimates and diagnostics, with columns:
#'     \emph{Estimate}, \emph{Std. Error} (pooled), 95\% CI bounds
#'     (\emph{2.5\%}, \emph{97.5\%}), degrees of freedom (\emph{DF}), and fraction
#'     of missing information (\emph{FMI}).}
#'   \item{m}{Number of draws (imputations) used in pooling.}
#'
#' @details Estimates are pooled using Rubin's rules: the pooled estimate is the mean
#' of draw-specific estimates, and the pooled variance is the sum of the within-imputation
#' variance and the scaled between-imputation variance. 95\% confidence intervals use a
#' t-distribution with the pooled degrees of freedom. The fraction of missing information
#' is computed as \eqn{(1+1/M)\,B/T} for each parameter.
#'
#' @examples
#' \dontrun{
#' # Assume `mix_mod` is a fitted mixtureBayes survival model and `data` was used to fit it:
#' pool_res <- with.surv_mixtureBayes(mix_mod, data)
#' print(pool_res)
#' }
#'
#' @import survival
#' @export

with.surv_mixtureBayes <- function(object, data, formula = object$call$formula, dist = object$family, min_n = NULL, ...) {
  # Extract MCMC label samples (each row is a draw of z labels: 1 = true match)
  if (!is.null(object$m_samples)) {
    m_samples <- object$m_samples
  } else if (!is.null(object$m)) {
    m_samples <- object$m
  } else {
    stop("Mixture model object does not contain m_samples")
  }
  m_samples <- as.matrix(m_samples)
  m <- nrow(m_samples)
  if (m < 1) stop("No samples available in model object")
  # Use formula and distribution from object if not specified
  # if (is.null(formula)) formula <- object$call$formula

  if (missing(formula) || is.null(formula)) {
    if (!is.null(object$call$formula)) {
      formula <- eval(object$call$formula, envir = parent.frame())
    } else {
      stop("Formula not found: please provide it explicitly or ensure object$call$formula exists.")
    }
  }

  if (is.null(dist)) dist <- object$family

  # Early exit if Gamma is requested — survreg does not support it
  if (tolower(dist) == "gamma") {
    stop("Multiple imputation is not available for dist = 'gamma': survreg() does not support the gamma distribution.\nUse dist = 'weibull', 'lognormal', etc., or fit using Stan.")
  }

  # Compute default min_n = p + 1
  mf <- stats::model.frame(formula, data)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  p <- ncol(X)
  if (is.null(min_n)) {
    min_n <- p + 1
  }
  # Prepare storage for coefficients and variances
  coefs_list <- list()
  vars_list <- list()
  used <- 0
  for (i in seq_len(m)) {
    z <- m_samples[i, ]
    if (length(z) != nrow(data)) {
      stop("Length of z does not match number of rows in data")
    }
    idx <- which(z == 1)
    if (length(idx) < min_n) next  # skip if too few matches
    data_sub <- data[idx, , drop = FALSE]
    # fit <- try(survival::survreg(formula = Surv(y, status) ~  X2, data = data_sub, dist = dist, ...), silent = TRUE)
    fit <- try(survival::survreg(formula = formula, data = data_sub, dist = dist, ...), silent = TRUE)

    if (inherits(fit, "try-error")) next
    beta <- stats::coef(fit)
    vc <- try(stats::vcov(fit), silent = TRUE)
    if (inherits(vc, "try-error") || is.null(vc)) next
    if (any(is.na(beta))) next
    used <- used + 1
    coefs_list[[used]] <- beta
    vars_list[[used]] <- diag(vc)[names(beta)]
  }

  if (used == 0) stop("No successful fits were obtained (try lowering min_n)")
  # Combine into matrices (rows = imputations, cols = parameters)
  coef_mat <- do.call(rbind, coefs_list)
  var_mat <- do.call(rbind, vars_list)
  M_eff <- nrow(coef_mat)
  # Pooled estimate (mean) and variances
  qbar <- colMeans(coef_mat)
  W <- colMeans(var_mat)
  B <- apply(coef_mat, 2, stats::var)
  Tvar <- W + (1 + 1/M_eff) * B
  std_error <- sqrt(Tvar)
  # Degrees of freedom (Rubin 1987 formula)
  df <- (M_eff - 1) * (1 + W/B)^2
  df[!is.finite(df)] <- Inf
  # 95% t-based confidence intervals
  alpha <- 0.05
  tval <- stats::qt(1 - alpha/2, df)
  ci_lower <- qbar - tval * sqrt(Tvar)
  ci_upper <- qbar + tval * sqrt(Tvar)
  # Fraction of missing information
  fmi <- ((1 + 1/M_eff) * B) / Tvar
  # Assemble output data frame
  out_df <- data.frame(
    Estimate = qbar,
    "Std. Error" = std_error,
    "2.5 %" = ci_lower,
    "97.5 %" = ci_upper,
    DF = df,
    FMI = fmi,
    check.names = FALSE
  )
  rownames(out_df) <- names(qbar)
  result <- list(call = object$call, pooled = out_df, m = M_eff)
  class(result) <- "mi_link_pool"
  return(result)
}

#' Print method for mi_link_pool objects
#'
#' @description
#' Print pooled coefficient estimates from multiple-imputation pooling
#' based on posterior linkage draws.
#'
#' @param x An object of class \code{mi_link_pool}.
#' @param ... Additional arguments passed to \code{stats::printCoefmat()}.
#'
#' @returns Invisibly returns \code{x}.
#'
#' @export
print.mi_link_pool <- function(x, ...) {
  cat("Pooled survival mixture model coefficients:\n")
  mat <- as.matrix(x$pooled)
  stats::printCoefmat(mat, P.values = FALSE, has.Pvalue = FALSE, signif.stars = FALSE, ...)
}
