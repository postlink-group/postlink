#' Print Method for `adjMixBayes` Objects
#'
#' Provides a concise summary of the Bayesian adjustment object created by
#' \code{\link{adjMixBayes}}.
#'
#' @param x An object of class \code{adjMixBayes}.
#' @param digits Number of significant digits used for derived prior
#'   hyperparameters.
#' @param ... Additional arguments passed to methods.
#' @return Invisibly returns \code{x}.
#'
#' @details
#' This method inspects the reference-based data environment to report the number
#' of linked records without copying the full dataset. It safely handles cases
#' where the linked data is unspecified (NULL). It also prints the user-specified
#' priors (marking entries that the fitting functions do not recognise, such
#' as \code{intercept} without its component number, as ignored), outlines the
#' defaults that will be used, and describes the match
#' probability model implied by \code{m.formula}, \code{m.rate} and
#' \code{safe.matches}, applying the same precedence rules as the fitting
#' functions (an explicit \code{theta} or \code{gamma_intercept} prior takes
#' precedence over \code{m.rate}; with covariates in \code{m.formula}, a
#' \code{gamma_intercept} prior also over a \code{theta} prior, and
#' \code{gamma_slope} over its alias \code{gamma}). When a shape parameter of
#' the beta distribution implied by \code{m.rate} and \code{m.rate.sd} is below
#' 1, a note gives the bound on \code{m.rate.sd} below which both exceed 1,
#' rounded down (see \code{m.rate.sd} in \code{\link{adjMixBayes}}). With covariates in \code{m.formula}, the
#' intercept prior and the prior on \eqn{\theta}{theta} implied by
#' \code{m.rate} refer to a record with average linkage covariates: the fit
#' centres them at their mean over the analysed records not flagged as safe
#' matches (over all analysed records when every record is a safe match) and
#' stores this mean as \code{z_center} in the fit. When the object holds the
#' linked data, the corresponding mean over the records of \code{linked.data}
#' is shown; it can differ slightly from \code{z_center} when records are
#' dropped at fit time (e.g. by \code{subset} or because of missing outcomes).
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
#' adjMixBayes(lifem_small, m.formula = ~ commf + comml, m.rate = 0.05,
#'             m.rate.sd = 0.02, safe.matches = hndlnk)
#'
#' @export
print.adjMixBayes <- function(x, digits = 3, ...) {
 cat("\n Adjustment Object: Bayesian Mixture \n")

 NextMethod("print")

 cat("\n* Priors:\n")
 if (!is.null(x$priors) && length(x$priors) > 0) {
  for (p in names(x$priors)) {
   # entries that the fitting functions do not recognise (e.g. "intercept"
   # without its component number) are ignored
   unused <- if (is.null(.key_menu(p))) "   [not recognised: ignored]" else ""
   cat(sprintf("      %-10s : %s%s\n", p, format(x$priors[[p]]), unused))
  }
  cat("    (Unspecified parameters will use defaults below)\n")
 } else {
  cat("\n    Status:       None specified; the defaults below are used.\n")
 }

 fmt <- function(v) format(v, digits = digits)
 use_logistic <- !.is_intercept_only_formula(x$m.formula)
 user <- if (is.list(x$priors)) names(x$priors) else character(0)
 user_theta <- "theta" %in% user
 user_gi <- use_logistic && "gamma_intercept" %in% user
 sd_val <- if (!is.null(x$m.rate.sd)) x$m.rate.sd else 0.1
 mrate_used <- !is.null(x$m.rate) && !user_theta && !user_gi
 bp <- if (mrate_used) tryCatch(mrate_to_beta(x$m.rate, sd_val), error = function(e) NULL) else NULL
 # records over which the fit centres the linkage covariates: those not
 # flagged as safe matches, or all records when every record is safe
 zc <- if (use_logistic) .adj_z_center(x) else NULL
 all_safe <- if (!is.null(zc)) isTRUE(attr(zc, "all_safe")) else
  length(x$safe.matches) > 0L && isTRUE(all(x$safe.matches == 1))

 cat("\n")
 cat("    Defaults applied during fitting (for any unspecified):\n")
 cat("      Intercept:   intercept1, intercept2 ~ normal(0,10)\n")
 cat("      GLM Slopes:  beta1, beta2 ~ normal(0,5) [binomial: beta1 ~ normal(0,2.5)]\n")
 cat("      Surv Slopes: beta1, beta2 ~ normal(0,5) [weibull: normal(0,2)]\n")
 cat("      Dispersion:  gaussian sigma1, sigma2 ~ cauchy(0,2.5); gamma phi1, phi2 ~ gamma(2,0.1)\n")
 cat("      Surv Shape:  gamma phi1, phi2 ~ exponential(1); weibull shape1, shape2 ~ gamma(2,1),\n")
 cat("                   weibull scale1, scale2 ~ gamma(2,1)\n")
 if (use_logistic) {
  int_txt <- if (user_gi) {
   if (user_theta) "user-specified (the theta prior is not used)" else "user-specified"
  } else if (user_theta) {
   tb <- tryCatch(.parse_prior("theta", x$priors$theta, .prior_menu$theta)$args, error = function(e) NULL)
   if (is.null(tb)) "derived from the theta prior" else {
    mom <- logit_beta_moments(tb[1], tb[2])
    sprintf("normal(%s,%s) [logit moments of the theta prior]", fmt(mom[["mu"]]), fmt(mom[["sd"]]))
   }
  } else if (!is.null(bp)) {
   mom <- logit_beta_moments(bp$alpha, bp$beta)
   sprintf("normal(%s,%s) [from m.rate]", fmt(mom[["mu"]]), fmt(mom[["sd"]]))
  } else {
   "normal(0,1.81) [logit of beta(1,1)]"
  }
  slope_txt <- if (all(c("gamma_slope", "gamma") %in% user)) {
   "user-specified (gamma_slope; its alias gamma is not used)"
  } else if (any(c("gamma_slope", "gamma") %in% user)) {
   "user-specified"
  } else {
   "normal(0,2.5)"
  }
  cat("      Match Model: gamma intercept ~ ", int_txt, "; slopes ~ ", slope_txt, "\n", sep = "")
  cat("                   (the intercept refers to a record with average linkage covariates: they\n", sep = "")
  if (all_safe) {
   cat("                   are centred at their mean over all records, as every record is flagged\n",
       "                   as a safe match)\n", sep = "")
  } else {
   cat("                   are centred at their mean over the records not flagged as safe matches)\n", sep = "")
  }
 } else if (user_theta) {
  cat("      Mix Weight:  theta: user-specified (see above)\n")
 } else if (!is.null(bp)) {
  cat("      Mix Weight:  theta ~ beta(", fmt(bp$alpha), ",", fmt(bp$beta), ") [from m.rate]\n", sep = "")
 } else {
  cat("      Mix Weight:  theta ~ beta(1,1)\n")
 }

 cat("\n* Match Probability Model:")
 if (use_logistic) {
  f_text <- deparse(x$m.formula)
  if (length(f_text) > 1) f_text <- paste(f_text[1], "...")
  cat("\n    Formula:              ", f_text, " (logistic regression)", sep = "")
 } else {
  cat("\n    Formula:               ~1 (constant match probability)")
 }

 if (!is.null(x$m.rate)) {
  if (mrate_used) {
   cat("\n    Prior Mismatch Rate:  ", fmt(x$m.rate), " (prior SD of theta ", fmt(sd_val), ")", sep = "")
   if (!is.null(bp)) {
    if (use_logistic) {
     mom <- logit_beta_moments(bp$alpha, bp$beta)
     q <- stats::plogis(mom[["mu"]] + c(0, -1, 1) * stats::qnorm(0.975) * mom[["sd"]])
    } else {
     q <- stats::qbeta(c(0.5, 0.025, 0.975), bp$alpha, bp$beta)
    }
    cat("\n    Implied Prior on theta: median ", fmt(q[1]), ", 95% interval (", fmt(q[2]), ", ",
        fmt(q[3]), ")", sep = "")
    if (use_logistic) {
     over <- if (all_safe) {
      "all records\n      (every record is flagged as a safe match)"
     } else {
      "the records not flagged as safe matches"
     }
     # the mean over linked.data; the fit averages over the analysed records
     where <- if (is.null(zc)) {
      "(computed when the model is fitted and stored as fit$z_center: no linked data are stored)"
     } else {
      paste0("(in linked.data: ", paste0(names(zc), " = ", vapply(zc, fmt, ""), collapse = ", "),
             "; fit$z_center holds the mean over the analysed records)")
     }
     cat("\n      for a record with average linkage covariates over ", over,
         "\n      ", where, sep = "")
    }
    # bound on m.rate.sd rounded down (see .mrate_dispersed_message())
    low <- .shape_below_one(c(bp$alpha, bp$beta))
    if (any(low)) {
     bound <- .mrate_sd_bound_text(x$m.rate)
     if (use_logistic) {
      psd <- .logitnormal_moments(mom[["mu"]], mom[["sd"]])[["sd"]]
      cat("\n    Note: on the probability scale this prior has standard deviation ", fmt(psd),
          " (rather than\n          m.rate.sd = ", fmt(sd_val), "); with m.rate.sd = ", bound, " it would be ",
          fmt(.mrate_logit_sd(x$m.rate, as.numeric(bound))), ".", sep = "")
     } else {
      where <- if (all(low)) "0 and 1" else if (low[2]) "1" else "0"
      cat("\n    Note: a shape parameter of this beta prior is below 1, so its density is unbounded",
          "\n          at ", where, "; with m.rate.sd at most ", bound, " both shape parameters exceed 1.",
          sep = "")
     }
    }
   }
  } else {
   cat("\n    Prior Mismatch Rate:  ", fmt(x$m.rate), " (not used: an explicit ",
       if (user_gi) "gamma_intercept" else "theta", " prior was supplied)", sep = "")
  }
 } else {
  cat("\n    Prior Mismatch Rate:   None specified (Estimated from data)")
 }

 if (!is.null(x$safe.matches)) {
  n_safe <- sum(x$safe.matches, na.rm = TRUE)
  n_all <- length(x$safe.matches)
  pct_safe <- if (n_all > 0) (n_safe / n_all) * 100 else 0
  cat("\n    Safe Matches:         ", format(n_safe, big.mark = ","),
      sprintf(" (%.1f%%)", pct_safe), sep = "")
 } else {
  cat("\n    Safe Matches:          None specified")
 }

 cat("\n\n")
 invisible(x)
}
