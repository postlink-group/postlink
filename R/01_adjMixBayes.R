#' Secondary Analysis Constructor Based on Bayesian Mixture Modeling
#'
#' Specifies the linked data and information on the underlying record linkage
#' process for a mismatch error adjustment using a Bayesian framework based on
#' mixture modeling as developed by Gutman et al. (2016). This framework uses
#' a mixture model for pairs of linked records whose two components reflect distributions
#' conditional on match status, i.e., correct match or mismatch.
#' Posterior inference is carried out via data augmentation or multiple
#' imputation.
#'
#' @param linked.data A data.frame containing the linked dataset. It is required
#'   for fitting models with \code{\link{plglm}()} or \code{\link{plsurvreg}()},
#'   which take the outcome, the covariates, the linkage covariates of
#'   \code{m.formula} and the safe matches from it.
#' @param priors A named \code{list} (or \code{NULL}) of prior specifications
#'   written as distribution strings (e.g. \code{"normal(0, 5)"}). These strings
#'   are checked when the object is created and parsed into numeric
#'   hyperparameters when the model is fitted. Any missing entries are
#'   automatically filled with defaults during the model fitting phase; see the
#'   section \emph{Prior distributions} of \code{\link{glmMixBayes}} and
#'   \code{\link{survregMixBayes}} for the supported distributions and the
#'   defaults. Priors given at fit time, e.g.
#'   \code{plglm(..., priors = list(...))} or \code{control = list(priors = ...)},
#'   replace the ones stored in the adjustment object.
#'
#'   For all models, intercept and slope priors are decoupled. Use
#'   \code{intercept1}/\code{intercept2} for the intercept of each component and
#'   \code{beta1}/\code{beta2} for the slope coefficients. This is particularly
#'   useful when the mismatch component slopes should be shrunk toward zero while
#'   allowing the intercept to remain unrestricted. Component-specific entries
#'   need the component number: names such as \code{intercept} or \code{beta}
#'   are not recognised and are ignored with a warning, given when the object
#'   is created (or, for entries given at fit time or added to the object
#'   later, when the model is fitted). When the design matrix of the outcome model has no intercept
#'   column (e.g. \code{y ~ 0 + x}), every column receives the
#'   \code{beta1}/\code{beta2} prior and \code{intercept1}/\code{intercept2}
#'   are not used (with a warning at fit time).
#'
#' @param m.formula A one-sided formula specifying the model for the conditional
#'   probability of a correct match given auxiliary linkage variables found in
#'   \code{linked.data} (e.g., \code{~ commf + comml}). The default (\code{~1})
#'   assumes a constant match probability \eqn{\theta}{theta} across all
#'   observations. When covariates are supplied, the match probability of each
#'   observation is modeled via a logistic regression
#'   \eqn{\theta_n = \mathrm{inv\_logit}(Z_n \gamma)}{theta_n = inv_logit(Z_n * gamma)}.
#'   The formula must keep its intercept and must not contain \code{offset()}
#'   terms. The \eqn{\gamma}{gamma} vector is decoupled into an intercept, with
#'   prior \code{gamma_intercept} (by default derived from \code{m.rate} or from
#'   the \code{theta} prior), and slopes, with prior \code{gamma_slope} (default
#'   \code{normal(0, 2.5)}); both must be \code{normal(mu, sd)}. The linkage
#'   covariates are centred at fit time at their mean over the records not
#'   flagged as safe matches, so that the \code{gamma_intercept} prior, whether
#'   derived or given, refers to a record with average linkage covariates (the
#'   centre is stored as \code{z_center} in the fit, and the reported
#'   \eqn{\gamma}{gamma} refers to the original covariates). For the default
#'   slope prior to be sensible, continuous covariates should be approximately on
#'   a unit scale; otherwise standardize them or supply
#'   \code{priors = list(gamma_slope = "normal(0, sd)")}. Records whose linkage
#'   covariates are missing are dropped at fit time with a warning.
#'
#' @param m.rate An optional numeric scalar strictly between 0 and 1 giving a
#'   prior estimate of the proportion of false links (mismatch rate) among the
#'   records that are not flagged in \code{safe.matches}.
#'   Together with \code{m.rate.sd}, it informs the prior on the match
#'   probability \eqn{\theta}{theta} by moment matching on the probability
#'   scale: \eqn{E[\theta] = 1 - \mathrm{m.rate}}{E[theta] = 1 - m.rate} and
#'   \eqn{SD[\theta] = \mathrm{m.rate.sd}}{SD[theta] = m.rate.sd}; the
#'   mismatch rate is then estimated from the data under this prior. These
#'   moments hold exactly with the default \code{m.formula = ~1}; with
#'   covariates in \code{m.formula} the intercept of \eqn{\gamma}{gamma}
#'   receives a normal prior on the logit scale (see \code{m.rate.sd}), whose
#'   probability-scale mean and standard deviation are close to these values
#'   only for a moderate \code{m.rate.sd} (for \code{m.rate = 0.05} the
#'   probability-scale standard deviation is about 0.077 for
#'   \code{m.rate.sd = 0.0475} and 0.27 for the default 0.1). This
#'   differs from \code{\link{adjMixture}()}, where \code{m.rate} is a
#'   constraint: \code{plglm()} and \code{plcoxph()} use it as an upper bound
#'   on the estimated mismatch rate (with covariates in \code{m.formula}, at
#'   the mean of the linkage covariates over the records not flagged as safe
#'   matches), and \code{plctable()} uses it as the known mismatch rate. Here,
#'   with covariates in \code{m.formula}, the prior derived from \code{m.rate}
#'   applies to the match probability of a record with average linkage
#'   covariates (the mean over the records not flagged as safe matches, the
#'   point used by \code{adjMixture()}): it is placed on the intercept of
#'   \eqn{\gamma}{gamma} for the centred covariates. A \code{theta} entry (or, with covariates, a
#'   \code{gamma_intercept} entry) in \code{priors} takes precedence over
#'   \code{m.rate}, with a message; with covariates, a \code{gamma_intercept}
#'   entry also takes precedence over a \code{theta} entry (with a message).
#'
#' @param m.rate.sd A positive numeric value (default \code{0.1}) giving the
#'   prior standard deviation of \eqn{\theta}{theta} on the probability scale
#'   that accompanies \code{m.rate}; when \code{m.rate} is used (no explicit
#'   \code{theta} or, with covariates in \code{m.formula},
#'   \code{gamma_intercept} prior), it must satisfy
#'   \eqn{\mathrm{m.rate.sd}^2 < \mathrm{m.rate}(1-\mathrm{m.rate})}{m.rate.sd^2 < m.rate * (1 - m.rate)}.
#'   With the default \code{m.formula = ~1}, \code{(m.rate, m.rate.sd)} is
#'   converted to a \code{beta(a, b)} prior on \eqn{\theta}{theta}; with
#'   covariates in \code{m.formula}, the intercept of \eqn{\gamma}{gamma}
#'   receives a normal prior with the exact logit-scale mean and standard
#'   deviation of that beta distribution. Smaller values give a more
#'   concentrated prior. For small mismatch rates the default is large: a
#'   shape parameter of the beta distribution is then below 1, so that its
#'   density is unbounded at 0 or 1 (with covariates, the normal prior on the
#'   logit scale is then much more dispersed than \code{m.rate.sd} suggests:
#'   for \code{m.rate = 0.05} and the default \code{m.rate.sd}, its
#'   probability-scale standard deviation is about 0.27). Both shape
#'   parameters exceed 1, and the density is bounded, when \code{m.rate.sd} is
#'   below
#'   \eqn{\sqrt{\min(m^2(1-m)/(1+m),\, m(1-m)^2/(2-m))}}{sqrt(min(m^2 (1 - m) / (1 + m), m (1 - m)^2 / (2 - m)))}
#'   with \eqn{m}{m} = \code{m.rate} (for \code{m.rate = 0.05}, at most
#'   0.04755); a warning is issued at fit time otherwise, which states this
#'   bound rounded down. A bounded density does not mean a mode
#'   near \code{1 - m.rate}: near the bound the density of the beta prior for
#'   a small \code{m.rate} rises up to \eqn{\theta = 1}{theta = 1} or very
#'   close to it (for \code{m.rate = 0.05} the prior is then about
#'   \code{beta(19, 1)}), and a mode near
#'   \code{1 - m.rate} needs a smaller value (e.g. \code{m.rate.sd = 0.02} for
#'   \code{m.rate = 0.05}, which gives a mode of about 0.96).
#'   \code{m.rate.sd} has no effect without \code{m.rate}.
#'
#' @param safe.matches An optional logical (or 0/1) vector of length
#'   \code{nrow(linked.data)} indicating records that are known to be correct
#'   matches (\code{TRUE}), for example hand-linked records. These records are
#'   always assigned to the correct-match component and do not inform the
#'   match probability model. Alternatively, the name of such a column in
#'   \code{linked.data} can be supplied, unquoted or as a character string.
#'
#' @return An object of class \code{c("adjMixBayes", "adjustment")}. To minimize
#' memory overhead, the underlying \code{linked.data} is stored by reference
#' within an environment inside this object.
#'
#' @details
#' \code{linked.data} must be supplied for the adjustment object to be used
#' with \code{plglm()} or \code{plsurvreg()}: the fitting methods take the
#' model variables from it and align the \code{m.formula} covariates and
#' \code{safe.matches} with the rows of the outcome model through its row
#' names; they refuse to fit an object created without data. Weights and
#' offsets are not supported by the Bayesian mixture models: \code{offset()}
#' terms in the model formula of \code{plglm()} or \code{plsurvreg()}, and
#' \code{weights} or \code{offset} arguments, raise an error.
#'
#' The fitted objects use the names of the frequentist mixture fits of
#' \code{\link{adjMixture}()}: \code{coefficients} for the outcome model of the
#' correct matches, \code{m.coefficients} for the mismatch-indicator model of
#' \code{m.formula} (on the logit scale of the mismatch probability, with the
#' same sign as in \code{adjMixture()}) and \code{match.prob} for the
#' probability that each record is a correct match; the outcome model of the
#' mismatches (component 2) is stored as \code{coefficients2}. The section
#' \emph{Correspondence with adjMixture()} of \code{\link{glmMixBayes}} and
#' \code{\link{survregMixBayes}} lists all correspondences and differences
#' (e.g. the meaning of \code{m.rate}).
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
#' # Construct the Bayesian mixture adjustment object. The priors are on the
#' # scale of the outcome analysed with lifem (age at death in years, regressed
#' # on a polynomial of the year of birth, as in ?plglm): the default priors
#' # suit outcomes of order one. beta2 = "normal(0, 0.01)" says that the
#' # mismatches show no trend in the year of birth.
#' age_priors <- list(
#'   intercept1 = "normal(60, 20)", intercept2 = "normal(60, 20)",
#'   beta1 = "normal(0, 100)", beta2 = "normal(0, 0.01)"
#' )
#' adj_bayes2 <- adjMixBayes(
#'   linked.data = lifem_small,
#'   priors = c(age_priors, theta = "beta(2, 2)")
#' )
#'
#' class(adj_bayes2)
#'
#' # Known correct matches, linkage covariates for the match probability and a
#' # prior estimate of the mismatch rate (with a prior SD small enough for a
#' # 5% mismatch rate)
#' adj_bayes3 <- adjMixBayes(
#'   linked.data = lifem_small,
#'   priors = age_priors,
#'   m.formula = ~ commf + comml,
#'   m.rate = 0.05,
#'   m.rate.sd = 0.02,
#'   safe.matches = hndlnk
#' )
#'
#' adj_bayes3
#'
#' @seealso
#' * [plglm()] for generalized linear regression modeling
#' * [plsurvreg()] for parametric survival modeling
#'
#' @references
#' Gutman, R., Sammartino, C., Green, T., & Montague, B. (2016). Error adjustments for file
#' linking methods using encrypted unique client identifier (eUCI) with application to recently
#' released prisoners who are HIV+. \emph{Statistics in Medicine}, 35(1), 115–129. \doi{10.1002/sim.6586}
#'
#' @export
adjMixBayes <- function(linked.data = NULL, priors = NULL,
                        m.formula = ~1, m.rate = NULL, m.rate.sd = 0.1,
                        safe.matches = NULL) {
 # 1. Validate linked.data
 if (!is.null(linked.data)) {
  if (is.environment(linked.data) || is.list(linked.data)) {
   linked.data <- tryCatch(as.data.frame(linked.data), error = function(e) {
    stop("'linked.data' must be a data.frame or coercible to one.", call. = FALSE)
   })
  } else if (!is.data.frame(linked.data)) {
   stop("'linked.data' must be a data.frame, list, or environment.", call. = FALSE)
  }
 }

 # 2. Validate priors: a named list of prior strings; recognised entries must
 #    use a supported distribution with valid arguments. The entries reported
 #    as not recognised are recorded, so that the fit ignores them without
 #    repeating the warning (entries added to the object later are reported
 #    at fit time)
 unknown <- .check_prior_list(priors)
 if (length(unknown) > 0L) attr(priors, "unrecognised") <- unknown

 # 3. Validate m.formula (one-sided, with an intercept and no offset;
 #    covariates must be in linked.data)
 if (is.null(m.formula)) m.formula <- ~1
 if (!inherits(m.formula, "formula")) {
  stop("'m.formula' must be a formula object.", call. = FALSE)
 }
 if (length(m.formula) != 2L) {
  stop("'m.formula' must be a one-sided formula (e.g., ~ z1 + z2).", call. = FALSE)
 }
 m_vars <- all.vars(m.formula)
 if ("." %in% m_vars) stop("Usage of '.' in 'm.formula' is not supported.", call. = FALSE)
 m_terms <- stats::terms(m.formula)
 if (attr(m_terms, "intercept") == 0L) {
  stop("'m.formula' must include an intercept (the baseline match probability); ",
       "remove '- 1' / '+ 0'.", call. = FALSE)
 }
 if (!is.null(attr(m_terms, "offset"))) {
  stop("offset() terms are not supported in 'm.formula'.", call. = FALSE)
 }
 if (!is.null(linked.data) && !.is_intercept_only_formula(m.formula)) {
  missing_vars <- setdiff(m_vars, names(linked.data))
  if (length(missing_vars) > 0L) {
   stop("The following variable(s) in 'm.formula' are not found in 'linked.data': ",
        paste(missing_vars, collapse = ", "), call. = FALSE)
  }
 }

 # 4. Validate m.rate and m.rate.sd (probability-scale SD of theta); whether
 #    the pair defines a beta prior matters only when m.rate is used, i.e.
 #    without an explicit theta (or, with covariates in m.formula,
 #    gamma_intercept) prior
 if (is.null(m.rate.sd)) {
  stop("'m.rate.sd' must be a single positive numeric value.", call. = FALSE)
 }
 explicit <- c("theta", if (!.is_intercept_only_formula(m.formula)) "gamma_intercept")
 mrate_used <- !any(vapply(explicit, function(k) is.list(priors) && !is.null(priors[[k]]), logical(1)))
 .validate_mrate(m.rate, m.rate.sd, check_pair = mrate_used)
 if (is.null(m.rate) && !missing(m.rate.sd)) {
  warning("'m.rate.sd' has no effect unless 'm.rate' is supplied.", call. = FALSE)
 }

 # 5. Resolve safe.matches (NSE within linked.data, then standard evaluation)
 safe_matches_eval <- NULL
 safe_expr <- substitute(safe.matches)
 if (!is.null(safe_expr)) {
  safe_matches_eval <- tryCatch(
   eval(safe_expr, linked.data, enclos = parent.frame()),
   error = function(e) NULL
  )
  if (is.null(safe_matches_eval)) {
   safe_matches_eval <- tryCatch(
    eval(safe_expr, envir = parent.frame()),
    error = function(e) {
     stop("Could not find object '", deparse(safe_expr),
          "' in 'linked.data' or the environment.", call. = FALSE)
    }
   )
  }
 }
 if (!is.null(safe_matches_eval)) {
  # the name of a column, given as a character string
  if (is.character(safe_matches_eval) && length(safe_matches_eval) == 1L &&
      !is.null(linked.data) && safe_matches_eval %in% names(linked.data)) {
   safe_matches_eval <- linked.data[[safe_matches_eval]]
  }
  if (is.numeric(safe_matches_eval) && all(safe_matches_eval %in% c(0, 1))) {
   safe_matches_eval <- safe_matches_eval == 1
  }
  if (!is.logical(safe_matches_eval)) {
   stop("'safe.matches' must be a logical (or 0/1) vector of length nrow(linked.data), ",
        "or the name of such a column of 'linked.data'.", call. = FALSE)
  }
  if (anyNA(safe_matches_eval)) {
   stop("'safe.matches' must not contain NA values.", call. = FALSE)
  }
  if (!is.null(linked.data) && length(safe_matches_eval) != nrow(linked.data)) {
   stop("'safe.matches' must have the same length as 'nrow(linked.data)' (",
        nrow(linked.data), ").", call. = FALSE)
  }
  safe_matches_eval <- as.vector(safe_matches_eval)
 }

 # 6. Construct and return the S3 object with reference semantics
 data_ref <- new.env(parent = emptyenv())
 data_ref$data <- linked.data

 out <- structure(
  list(
   data_ref     = data_ref,
   priors       = priors,
   m.formula    = m.formula,
   m.rate       = m.rate,
   m.rate.sd    = m.rate.sd,
   safe.matches = safe_matches_eval
  ),
  class = c("adjMixBayes", "adjustment")
 )

 return(out)
}
