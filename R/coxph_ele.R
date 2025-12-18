#' Fit a CoxPH Model Assuming Exchangeable Linkage Errors
#' @import survival
#' @import nleqslv
#'
#' @description
#' Fit Cox proportional hazards regression adjusted for mismatched data based
#' on the approach developed in Vo et al., 2024 assuming exchangeable linkage
#' errors. Block-wise mismatch rates are assumed to be known.
#'
#' @param formula a formula object for the outcome model, the response should be
#' provided using the `Surv` function and the covariates should be separated by + signs.
#' @param data a data.frame with linked data used in "formula" and "formula.m" (optional)
#' @param m.rate block-wise mismatch rates (should be a vector with length equal
#' to the number of blocks) - by default assume a single block.
#' @param blocks block indicators.
#' @param control an optional list variable to of control arguments including
#' "init.beta" for the initial outcome model coefficient estimates) - by
#' default is the naive estimator.
#' @param ... the option to directly pass "control" arguments
#'
#' @returns a list of results from the function called depending on the "family"
#' specified.
#' \item{coefficients}{the outcome model coefficient estimates}
#' \item{standard.errors}{the estimated standard errors}
#' \item{call}{the matched call}
#' \item{vcov}{the variance-covariance matrix}
#' \item{linear.predictors}{the linear predictors}
#'
#' @note
#' The references below discuss the implemented framework in more detail.
#'
#' @references Vo, T. H., Garès, V., Zhang, L. C., Happe, A., Oger, E.,
#' Paquelet, S., & Chauvet, G. (2024). Cox regression with linked data.
#' Statistics in Medicine, 43(2), 296-314.\cr
#'
#' @export
coxph_ele <- function(formula, data,
                      m.rate, blocks,
                      control = list(init.beta = "default"),...){
 dcontrols <- list(...)
 if ("init.beta" %in% names(dcontrols)){
  init.beta <- dcontrols$init.beta
 } else {
  init.beta <- ifelse("init.beta" %in% names(control),
                      control$init.beta, "default")
 }

 if(missing(formula)){ stop("Error: a formula for the outcome
                            model is required")}
 if(!inherits(formula, "formula")){ stop("Error: formula should be
                                         a formula object")}

if(!(missing(data) | is.null(data)) && (!is.data.frame(data) & !is.list(data))){
  stop("Error: data should be a data.frame or list")}

 # Define Z, delta, and T based on formula and data (if provided)
 if (!missing(data)){
  mf <- model.frame(formula, data = data)
  delta <- mf[[1]][,2]
  Z <- mf[,-1]
  T <- mf[[1]][,"time"]
  if(attr(mf[[1]], "type") != "right"){
   stop("Error: censoring type other than right-censoring
        is not currently supported")
  }
  if(!is.Surv(mf[[1]])){
   stop(("Error: response should be a survival object"))}
 } else {
  mf <- model.frame(formula)
  delta <- mf[[1]][,2]
  Z <- mf[,-1]
  T <- mf[[1]][,"time"]
  if(attr(mf[[1]], "type") != "right"){
   stop("Error: censoring type other than right-censoring
        is not currently supported")
  }
  if(!is.Surv(mf[[1]])){
   stop(("Error: response should be a survival object"))}
 }
 n <- nrow(Z)
 p <- ncol(Z)

 if(any(is.na(Z)) | any(is.na(T))){"Error (formula):
  Cannot have a missing observations"}

 # Define initial beta for when we solve the AEE
 if (init.beta != "default"){
  init.beta <- as.vector(init.beta)
  if(length(init.beta) != p){
   warning("Default 'init.beta' used.
           'init.beta' should be a vector of length ", p)
   init.beta <- "default"}
 }
 if (init.beta == "default"){
  init.beta <- coxph(Surv(time = T, event = delta) ~ Z)$coefficients
 }

 # Define a_nu and blocks
 if(missing(blocks) | is.null(blocks)){
  blocks <- rep(1, n)
  warning("'blocks' argument is missing or NULL - assuming a single block for all observations.")
 }
 # TO-DO: handling NA values for observations
 #na_rows <- attr(mf, "na.action")
 #if(!is.null(na_rows) && length(blocks) == n + length(na_rows)){
 # blocks <- blocks[-na_rows]
 #}
 if(length(blocks) != n){
  stop(("Error: 'blocks' does not have a length equal to the number of
        observations."))
 }
 a_nu <- 1 - m.rate
 if(length(a_nu) == 1 && length(unique(blocks)) != 1){
  a_nu <- rep(a_nu, length(unique(blocks)))
 }

 # Define bias-adjusted estimating equation
 corr_fun <- function(mat){
  t1 <- sweep(mat, MARGIN = 1,  FUN = "*", STATS = (a_nu[blocks])^(-1))
  t2 <- as.matrix(aggregate(mat,by = list(blocks),
                            FUN = mean, drop = TRUE))[blocks,-1] *((a_nu[blocks])^(-1) - 1)
  return(t1 - t2)
 }

 hi_fun <- function(beta_est){
  Xstar <- corr_fun(Z)
  gstar <- corr_fun(exp(Z %*% beta_est))
  hstar <- corr_fun(sweep(Z, 1, exp(Z %*% beta_est), "*"))

  risk_set <- apply(as.matrix(T), 1, function(Ti) which(T >= Ti))
  gj <- unlist(lapply(risk_set, function(x) sum(gstar[x,])))
  hj <- do.call(rbind, lapply(risk_set, function(x) colSums(hstar[x,])))
  ratio <- hj/gj

  H_i <- (Xstar - ratio)[delta == 1,]
  return(H_i)
 }
 hbar_fun <- function(beta_est){(1/n)*colSums(hi_fun(beta_est))}

 # Obtain coefficient estimates
 res <- nleqslv(init.beta, hbar_fun, jacobian = TRUE)
 coef <- res$x

 # Obtain (sandwich-estimator) variance
 #nabla <- diag(res$jac)
 nabla <- res$jac
 #s2H <- (1/(n - 1)) * colSums((sweep(hi_fun(coef), 2, hbar_fun(coef)))^2)
 s2H <- (1/(n - 1)) * crossprod(sweep(hi_fun(coef), 2, hbar_fun(coef)))
 V2 <- s2H/n # assume a_nu is known for now
 #var <- nabla^(-1) %*% cbind(V2) %*% nabla^(-1) # variance per equations 8-9
 var <- solve(nabla) %*% V2 %*% t(solve(nabla))

 # Return results
 output <- list(coefficients = coef, standard.errors = diag(sqrt(var)),
                vcov = var, linear.predictors = X %*% coef, family = family)
 x <- append(output, match.call())
 names(x)[[length(x)]] <- "call"

 class(x) <- "coxph_ele"
 x
}
