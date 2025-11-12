#' Fit a GLM Assuming Exchangeable Linkage Errors
#' @import nleqslv
#' @description
#' Fit generalized linear regression models adjusted for mismatched data according
#' to the ELE approach (Chambers, 2009).
#' The function currently supports Gaussian, gamma, Poisson, and logistic (binary)
#' models. Block-wise mismatch rates are assumed to be known.
#'
#' @param formula a formula object for the outcome model, with the covariate(s) on
#' the right of "~" and the response on the left.
#' @param data a data.frame with linked data used in "formula" and "formula.m" (optional)
#' @param family the type of regression model ("gaussian" - default, "poisson",
#' "binomial", "gamma"). Standard link functions are used ("identity" for Gaussian,
#' "log" for Poisson and Gamma, and "logit" for binomial).
#' @param m.rate block-wise mismatch rates (should be a vector with length equal
#' to the number of blocks) - by default assume a single block.
#' @param blocks block indicators.
#' @param weight.matrix the type of weight matrix ("ratio-type",
#' "Lahiri-Larsen", "BLUE", or "all" (default))
#' @param control an optional list variable to of control arguments including
#' "init.beta" for the initial outcome model coefficient estimates) - by
#' default is the naive estimator when the weight matrix is ratio-type or
#' Lahiri-Larsen and is the Lahiri-Larsen estimator for the BLUE weight matrix.
#' @param ... the option to directly pass "control" arguments
#'
#' @returns a list of results from the function called depending on the "family" specified.
#' \item{coefficients}{the outcome model coefficient estimates}
#' \item{family}{the type of (outcome) regression model}
#' \item{standard.errors}{the estimated standard errors}
#' \item{call}{the matched call}
#' \item{dispersion}{the dispersion parameter estimate}
#' \item{vcov}{the variance-covariance matrix}
#' \item{linear.predictors}{the linear predictors}
#' \item{df.residual}{the residual degrees of freedom}
#'
#' @note
#' The references below discuss the implemented framework in more detail.
#'
#' @references Chambers, R. (2009). Regression analysis of probability-linked data.\cr
#'
#' @export
glm_ele <- function(formula, data = NULL, family = "gaussian",
                    m.rate, blocks, weight.matrix = "all",
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

 if(!missing(data) && (!is.data.frame(data) & !is.list(data))){
  stop("Error: data should be a data.frame or list")}

 if(!(family %in% c("gaussian", "poisson", "binomial", "gamma"))){
  stop("Error: the family should be gaussian, poisson, binomial, or gamma")}

 # Define X and y based on formula and data (if provided)
 if (!missing(data)){
  mf <- model.frame(formula, data=data)
  X <- model.matrix(formula, data=data)
  y <- model.response(model.frame(formula, data = data))
  if(family == "binomial"){
   y <- as.numeric(as.vector(y))
   if((sum(y > 1) != 0)) stop("Error: y values must be 0 or 1")
  } else {
   y <- as.numeric(y)
  }
  n <- nrow(X)
  p <- ncol(X)
 } else {
  mf <- model.frame(formula)
  X <- model.matrix(formula)
  y <- as.numeric(model.response(model.frame(formula)))
  if(family == "binomial" & (sum(y > 1) != 0)){
   stop("Error: currently only binary models are supported
        (# number of trials should be 1)")
  }
  n <- nrow(X)
  p <- ncol(X)
 }
 if(any(is.na(X)) | any(is.na(y))){"Error (formula):
  Cannot have a missing observations"}

 # Define initial beta for when we solve the AEE
 if (init.beta != "default"){
  init.beta <- as.vector(init.beta)
  if(length(init.beta) != p){
   warning("Default 'init.beta' used. 'init.beta' should
           be a vector of length ", p)
   init.beta <- "default"}
 }

 if (init.beta == "default"){
  if (family != "gamma"){
   init.beta <- coef(glm(y ~ X-1, family = family))
  } else{
   init.beta <- coef(glm(y ~ X-1, family = Gamma(link = "log")))
  }
 }

 # Define lambda, blocks, and weight.matrix (if "all")
 if(length(weight.matrix) == 1 && weight.matrix == "all"){
  weight.matrix = c("ratio-type", "Lahiri-Larsen", "BLUE")
 }
 if(missing(blocks)){
   blocks <- rep(1, n)
   warning("'blocks' argument is missing - assuming there is a single block.")
 }
 na_rows <- attr(mf, "na.action")
 if(!is.null(na_rows) && length(blocks) == n + length(na_rows)){
   blocks <- blocks[-na_rows]
 }
 if(length(blocks) != n){
  stop(("Error: 'blocks' does not have a length equal to the number of
        observations."))
 }
 lambda <- 1 - m.rate
 if(length(lambda) == 1 && length(unique(blocks)) != 1){
  lambda <- rep(lambda, length(unique(blocks)))
 }

 # Efficient matrix calculation (to avoid materializing Eq)
 IEq.M <- function(M, lambda, blocks, type = "Eq"){
  nq <- as.numeric(tapply(y, INDEX = blocks, FUN = length))
  alphaq <- 1-lambda
  if(type == "Eq"){
   c0q <- alphaq*(nq/(nq-1))
   c1q <- (1-alphaq - alphaq/(nq-1))
  }
  if(type == "Iq"){
   c0q <- -alphaq*(nq/(nq-1))
   c1q <- 1-(1-alphaq - alphaq/(nq-1))
  }
  if(!is.matrix(M) || ncol(M) == 1){
   prod <- M*c1q[blocks] + as.matrix(tapply(M, list(blocks), mean)[blocks]
                                     * c0q[blocks])
  } else{
   corr <- as.matrix(aggregate(M, by = list(blocks),
                               FUN = mean, drop = TRUE))[blocks,-1] * c0q[blocks]
   Mtilde <- sweep(M, MARGIN = 1,  FUN = "*", STATS = c1q[blocks])
   prod <- Mtilde + corr
  }
  return(prod)
 }

 # Estimate dispersion (sigma^2) for "gaussian" and "gamma" families
 sigma2 <- function(beta_hat, family){
  fq <- fq_fun(family, beta_hat)$fq
  num1 <- crossprod(y - fq)
  num2 <- -2*crossprod(fq, IEq.M(fq, lambda, blocks, type = "Iq"))
  num <- ifelse(num1+num2 < 0, num1/2, num1+num2)
  if(family == "gaussian"){
   sigma2 <- c(num/nrow(X))
  }
  if(family == "gamma"){
   sigma2 <- c(num/crossprod(fq))
  }
  return(sigma2)
 }

 # Estimate Var(Y_q^*) for BLUE weighting matrix and se calculation
 Var_y <- function(beta_hat, family){
  fq <- fq_fun(family, beta_hat)$fq
  fqbar <- as.matrix(tapply(fq, list(blocks), mean)[blocks])
  fq2bar <- as.matrix(tapply(fq^2, list(blocks), mean)[blocks])

  Vq <- (1-lambda[blocks])*(lambda[blocks]*(fq - fqbar)^2 + fq2bar - fqbar^2)

  if(family == "gaussian"){
   Exp <- sigma2(beta_hat, family)
  } else{
   if(family == "binomial"){
    v <- fq*(1-fq)
   }
   if(family == "poisson"){
    v <- fq
   }
   if(family == "gamma"){
    v <- sigma2(beta_hat, family)*fq^2
   }
   vsum <- as.matrix(tapply(v, list(blocks), sum)[blocks])
   nq <- as.numeric(tapply(y, INDEX = blocks, FUN = length))
   gamma <- (1-lambda)/(nq - 1)
   Exp <- (lambda[blocks] - gamma[blocks])*v + (gamma[blocks])*vsum
  }
  return(Exp + Vq)
 }

 # Define E(Y_q) = f_q(\beta)
 fq_fun <- function(family, beta_est){
  if(family == "gaussian"){
   fq <- X %*% beta_est
   dfq <- X
  }
  if(family == "binomial"){
   fq <- exp(X %*% beta_est)/(1+exp(X %*% beta_est))
   dfq <- sweep(X, MARGIN = 1,  FUN = "*", STATS = fq*(1-fq))
  }
  if(family %in% c("poisson", "gamma")){
   fq <- exp(X %*% beta_est)
   dfq <- sweep(X, MARGIN = 1,  FUN = "*", STATS = fq)
  }
  return(list(fq = fq, dfq = dfq))
 }

 # Define weighting matrix
 Gq_fun <- function(family, weight.matrix, beta_est){
  if(weight.matrix == "ratio-type"){
   Gq <- X
  }
  if(weight.matrix == "Lahiri-Larsen"){
   Gq <- IEq.M(X, lambda, blocks)
  }
  if(weight.matrix == "BLUE"){
   dfq <- fq_fun(family, beta_est)$dfq
   Gq <- sweep(IEq.M(dfq, lambda, blocks), MARGIN = 1,  FUN = "*",
               STAT = (Var_y(beta_est, family))^(-1))
  }
  return(Gq)
 }

 blocks <- match(blocks, sort(unique(blocks)))
 # Obtain results based on specified "weight.matrix"
 if("BLUE" %in% weight.matrix){
   weights.find <- unique(c(weight.matrix, "Lahiri-Larsen"))
 }
 all.weights <- c("ratio-type", "Lahiri-Larsen", "BLUE")
 weights.find <- all.weights[all.weights %in% weights.find]

 nwm <- length(weights.find)
 coef <- matrix(nr = nwm, nc = p, data = NA)
 var <- matrix(nr = nwm, nc = p, data = NA)
 covhat <- list()
 for (i in 1:nwm){
   wm <- weights.find[i]

   # Define bias-adjusted estimating equation
   AEE_fun <- function(beta_est){
     Gq <- Gq_fun(family, wm, beta_est)
     fq <- fq_fun(family, beta_est)$fq
     Eqfq <- IEq.M(fq, lambda, blocks)
     return(crossprod(Gq, y) - crossprod(Gq, Eqfq))
   }

   # Obtain coefficient estimates
   if(wm == "BLUE"){
     init.beta <- coef[which(weights.find == "Lahiri-Larsen"),]
   }
   res <- nleqslv(init.beta, AEE_fun, jacobian = TRUE)
   coef[i,] <- res$x

   # Obtain (sandwich-estimator) variance
   dH <- -res$jac
   Gq <- Gq_fun(family, wm, coef[i,])
   Vy <- Var_y(coef[i,], family)
   Vh <- crossprod(sweep(Gq, MAR = 1, FUN = "*", STAT = sqrt(Vy)))
   covhat[[i]] <- solve(dH) %*% Vh %*% t(solve(dH))
   var[i,] <- diag(covhat[[i]])
 }

 # Return results
 rownames(coef) <- weights.find
 rownames(var) <- weights.find
 colnames(coef) <- colnames(X)
 names(covhat) <- weights.find

 coef <- coef[weight.matrix,,drop=F]
 var <- var[weight.matrix,,drop=F]

 linear.predictors <- X %*% t(coef)
 colnames(linear.predictors) <- weight.matrix

 output <- list(coefficients = coef, standard.errors = sqrt(var),
                vcov = covhat, linear.predictors = linear.predictors,
                family = family)

 nwm <- length(weight.matrix)
 if(family %in% c("gamma", "gaussian")){
   output <- append(output, list(sapply(c(1:nwm),
                                        function(i) sigma2(coef[i,], family)),
                                 n - p - 1))
   names(output)[[(length(output)-1)]] <- "dispersion"
   names(output)[[(length(output))]] <- "df.residual"
   names(output$dispersion) <- weight.matrix
 }
 x <- append(output, match.call())
 names(x)[[length(x)]] <- "call"

 class(x) <- "glm_ele"
 x
}
