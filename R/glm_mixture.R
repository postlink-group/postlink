#' Fit a GLM Using a Mixture-Modeling Approach
#' @description
#' Fit generalized linear regression models adjusted for mismatched data.
#' The function currently supports Gaussian, gamma, Poisson, and logistic (binary)
#' models. Information about the underlying record linkage process can be
#' incorporated into the method if available (e.g., assumed overall mismatch rate,
#' safe matches, predictors of match status, or predicted probabilities of correct
#' matches).
#'
#' @param formula a formula object for the outcome model, with the covariate(s) on
#' the right of "~" and the response on the left.
#' @param data a data.frame with linked data used in "formula" and "formula.m" (optional)
#' @param family the type of regression model ("gaussian" - default, "poisson",
#' "binomial", "gamma"). Standard link functions are used ("identity" for Gaussian,
#' "log" for Poisson and Gamma, and "logit" for binomial).
#' @param m.formula a one-sided formula object for the mismatch indicator model, with the
#' covariates on the right of "~". The default is an intercept-only model corresponding
#' to a constant mismatch rate)
#' @param safe.matches an indicator variable for safe matches (TRUE : record can be treated as a
#' correct match and FALSE : record may be mismatched). The default is FALSE for all matches.
#' @param m.rate the assumed overall mismatch rate (a proportion between 0 and 1). If
#' not provided, no overall mismatch rate is assumed.
#' @param control an optional list variable to customize the initial parameter estimates
#' ("init.beta" for the outcome model and "init.gamma" for the mismatch indicator model),
#' estimated marginal density of the response ("fy"), maximum iterations for the
#' EM algorithm ("max.iter"), maximum iterations for the subroutine in the constrained
#' logistic regression function ("cmax.iter"), and convergence tolerance for
#' the termination of the EM algorithm ("tol").
#' @param ... the option to directly pass "control" arguments
#'
#' @returns a list of results from the function called depending on the "family" specified.
#' \item{coefficients}{the outcome model coefficient estimates}
#' \item{match.prob}{the posterior correct match probabilities for observations given parameter estimates}
#' \item{objective}{a variable that tracks the negative log pseudo-likelihood for all iterations of the EM algorithm.}
#' \item{family}{the type of (outcome) regression model}
#' \item{standard.errors}{the estimated standard errors}
#' \item{m.coefficients}{the correct match model coefficient estimates}
#' \item{call}{the matched call}
#' \item{wfit}{an internal-use object for the predict function}
#' \item{dispersion}{the dispersion parameter estimate}
#' \item{vcov}{the variance-covariance matrix}
#'
#' @examples
#' ## commonness score of first and last names used for linkage
#' mformula <- ~commf + comml
#' ## hand-linked records are considered "safe" matches
#' safematches <- ifelse(lifem$hndlnk =="Hand-Linked At Some Level", TRUE, FALSE)
#' ## overall mismatch rate in the data set is assumed to be ~ 0.05
#' mrate <- 0.05
#'
#' fit <- fit_mixture(age_at_death ~ poly(unit_yob, 3, raw = TRUE), data = lifem,
#'                    family = "gaussian", mformula, safematches, mrate)
#'
#' @note
#' The references below discuss the implemented framework in more detail. The standard
#' errors are estimated using the sandwich formula (Slawski et al., 2023).\cr
#'
#' @references Slawski, M., West, B. T., Bukke, P., Diao, G., Wang, Z., & Ben-David, E. (2023).
#' A General Framework for Regression with Mismatched Data Based on Mixture Modeling.
#' Under Review. < \doi{10.48550/arXiv.2306.00909} >\cr
#'
#' Slawski, M., Diao, G., Ben-David, E. (2021). A pseudo-likelihood approach to linear
#' regression with partially shuffled data. Journal of Computational and Graphical
#' Statistics. 30(4), 991-1003 < \doi{10.1080/10618600.2020.1870482} >
#'
#' @export
glm_mixture <- function(formula, data,
                        family = c("gaussian", "poisson", "binomial", "gamma"),
                        m.formula, safe.matches, m.rate,
                        control = list(init.beta = "default",
                                       init.gamma = "default",
                                       fy = "default",
                                       max.iter = 1000,
                                       tol = 1E-4, cmax.iter = 1000),...){
 # -----------------------------------------------------------------------------
 ### CONTROL ARGUMENTS ###
 # -----------------------------------------------------------------------------
 dcontrols <- list(...)
 if ("init.beta" %in% names(dcontrols)){
  init.beta <- dcontrols$init.beta
 } else {
  init.beta <- if("init.beta" %in% names(control)){control$init.beta} else{"default"}
 }

 if ("init.gamma" %in% names(dcontrols)){
  init.gamma <- dcontrols$init.gamma
 } else {
  init.gamma <- if("init.gamma" %in% names(control)){control$init.gamma} else{"default"}
 }

 if ("fy" %in% names(dcontrols)){
  fy <- dcontrols$fy
 } else {
  fy <- if("fy" %in% names(control)){control$fy} else{"default"}
 }

 if ("max.iter" %in% names(dcontrols)){
  max.iter <- dcontrols$max.iter
 } else {
  max.iter <- ifelse("max.iter" %in% names(control), control$max.iter, 1000)
 }

 if(max.iter < 2 | (floor(max.iter) != max.iter)){
  warning("Default max.iter used. max.iter should be an integer greater than 2")
  max.iter <- 1000}

 if ("cmax.iter" %in% names(dcontrols)){
  cmax.iter <- dcontrols$cmax.iter
 } else {
  cmax.iter <- ifelse("cmax.iter" %in% names(control), control$cmax.iter, 1000)
 }

 if(cmax.iter < 0 | (floor(cmax.iter) != cmax.iter)){
  warning("Default cmax.iter used. cmax.iter should be an integer")
  cmax.iter <- 1000}

 if ("tol" %in% names(dcontrols)){
  tol <- dcontrols$tol
 } else {
  tol <- ifelse("tol" %in% names(control), control$tol, 1E-4)
 }

 # -----------------------------------------------------------------------------
 ### MAIN ARGUMENTS ###
 # -----------------------------------------------------------------------------
 if(!missing(m.rate) && (m.rate <= 0 | m.rate >= 1)){
  stop("Error: assumed mismatch rate should be a proportion between 0 and 1")}

 if(missing(formula)){ stop("Error: a formula for the outcome model is required")}
 if(!inherits(formula, "formula")){ stop("Error: formula should be a formula object")}

if(!(missing(data) | is.null(data)) && (!is.data.frame(data) & !is.list(data))){
  stop("Error: data should be a data.frame or list")}

 family <- match.arg(family)
 if(!(family %in% c("gaussian", "poisson", "binomial", "gamma"))){
  stop("Error: the family should be gaussian, poisson, binomial, or gamma")}

 if(!missing(m.formula) && !inherits(m.formula, "formula")){
  stop("Error: m.formula should be a formula object")}

 if(!missing(safe.matches) && !missing(data)){
  val <- data[[deparse(substitute(safe.matches))]]
  if(!is.null(val)){
   safe.matches <- val
  }
 }

 if(!missing(safe.matches) && !is.logical(safe.matches)){
  stop("Error: safe.matches should be a logical object")}

 if (missing(data)){
  data <- NULL
 }

 # -----------------------------------------------------------------------------
 if (family == "gaussian"){
  x <- fit_mixture_gaussian(formula, data, family,
                            m.formula, safe.matches, m.rate,
                            init.beta, init.gamma, fy, max.iter, tol, cmax.iter)
 }

 if (family != "gaussian"){
  x <- fit_mixture_glm(formula, data, family,
                       m.formula, safe.matches, m.rate,
                       init.beta, init.gamma, fy, max.iter, tol, cmax.iter)
 }
 # -----------------------------------------------------------------------------
 x <- append(x, match.call())
 names(x)[[length(x)]] <- "call"

 class(x) <- "glm_mixture"
 x

}
