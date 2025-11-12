#' Spline least squares regression Based on Mixture Modeling
#' @import mgcv
#'
#' @description
#' This function performs least squares spline fitting accounting
#' for mismatched (x,y) pairs. Spline functionality is imported
#' from the {mgcv} package. Fitting is based on the Bayesian
#' formulation of the mixture model as described in Slawski et al. (2024).
#' Parameters are estimated using a vanilla mean-field variational approximation
#' of the posterior. The updates can be performed in closed form.
#' Currently only a constant mismatch rate model is supported.
#'
#' @param formula model formula of the form y ~ s(x), where x and y are
#' continuous variables contained in "data". This formula will be passed on
#' to calls of functions in the mgcv package, and can be used to specify
#' the number of knots.
#' @param data underlying data set
#' @param control an optional list variable with additional control arguments.
#  --- fy: estimated marginal density of the response
#  --- max.iter: maximum number of Variational Inference iterations (default: 100). Currently the algorithm is run for a fixed number of iterations).
#' @param ... other arguments to be passed to the gam() call in mgcv.
#'
#' @returns a list of results having the following components.
#' \item{B}{matrix of spline basis function evaluations from mgcv. Necessary for computing fitted values.}
#' \item{par_coef}{the parameters of the Normal distribution for the spline coefficients, i.e., its mean and covariance matrix.}
#' \item{par_sigmasq}{the parameters of the Inverse-Gamma distribution. (shape and rate)of the error variance sigmasq}
#' \item{par_tausq}{the parameters of the Inverse-Gamma distribution (shape and rate) of the variance of the (conditional) prior for the spline coefficients.}
#' \item{par_alpha}{the parameters of the Beta distribution for the mismatch rate alpha.}
#' \item{par_m}{the parameters of the Bernoulli distributions of the mismatch indicator (posterior probability of mismatch)}
#'
#' @export
spline_mixture <- function(formula, data,
                           control = list(max.iter = 100, fy = "default"), ...){
 if(missing(formula)){ stop("Error: a formula for the outcome model is required")}
 if(!inherits(formula, "formula")){ stop("Error: formula should be a formula object")}

 if(!missing(data) && (!is.data.frame(data) & !is.list(data))){
  stop("Error: data should be a data.frame or list")}
 if(missing(data)){data <- list()}

 # initialization with the naive estimator
 gam_naive <- gam(formula = formula, data = data, ...)
 gamobj <- gam(formula = formula, data = data, fit = FALSE)

 # extract quantities to be used later
 beta_naive <- coef(gam_naive)
 n <- length(gam_naive$fitted)
 # extract spline basis evaluations and smoothing matrix
 Sgam <- rbind(0, cbind(0, gamobj$smooth[[1]]$S[[1]]))
 Bgam <- as.matrix(gamobj$X)

 d <- ncol(Bgam)
 d_null <- sum(svd(Sgam)$d < sqrt(.Machine$double.eps))
 dstar <- d - d_null

 # obtain fy (currently only KDE)
 y <- gamobj$y

 fy <- control$fy
 if (!identical(fy, "default")){
  fy <- as.vector(fy)
  if(length(fy) != n){
   warning("Default 'fy' used. 'fy' should be a vector with a value
              for each observation")
   fy <- "default"}
 }
 if (identical(fy, "default")){
  kdeobj <- density(y)
  kdefun <- approxfun(kdeobj$x, kdeobj$y, method = "constant")
  fy <- kdefun(y)
 }

 # VI fitting
 VI_iter <- control$max.iter

 # initialize parameters
 a <- 1
 b <- 1
 par_beta <- list(mu = beta_naive, Prec = crossprod(Bgam))
 par_sigmasq <- list(shape = (n - d)/2, rate = sum(residuals(gam_naive)^2)/2)
 par_pi <- list(a = 1, b = 1) # prior parameters
 par_z <- numeric(n)
 # variance of the beta's
 par_tausq <- list(shape = dstar/2,
                   rate = sum(beta_naive * (Sgam %*% beta_naive)) / 2)
 #
 for(ii in 1:VI_iter){

  # [1] update z (Bernoulli distribution) given other variables

  E_log_pi_a <- digamma(par_pi$a) - digamma(sum(par_pi$a + par_pi$b))
  E_log_pi_b <- digamma(par_pi$b) - digamma(sum(par_pi$a + par_pi$b))

  log_prob_mm_unnorm <- E_log_pi_a + log(fy)


  E_reciproc_sigmasq <- par_sigmasq$shape / par_sigmasq$rate
  E_fit <- Bgam %*% par_beta$mu
  E_Quad <- colSums(t(Bgam) * solve(par_beta$Prec, t(Bgam))) + E_fit^2*E_reciproc_sigmasq
  E_log_reciproc_sigmasq <- (digamma(par_sigmasq$shape)) - log(par_sigmasq$rate)
  log_prob_cm_unnorm <- E_log_pi_b -0.5*(y^2 - 2*y * E_fit)*E_reciproc_sigmasq - 0.5*log(2*pi) + 0.5 * E_log_reciproc_sigmasq -0.5*E_Quad

  logpmax <- pmax(log_prob_mm_unnorm, log_prob_cm_unnorm)
  logsumexp <- logpmax + log(exp(log_prob_mm_unnorm - logpmax) +  exp(log_prob_cm_unnorm - logpmax))
  par_m <- exp(log_prob_mm_unnorm - logsumexp)


  # [2] update Beta distributions given z

  sum_par_m <- sum(par_m)
  par_pi$a <- a + sum_par_m
  par_pi$b <- b + n - sum_par_m

  # [3] update beta, sigma^2 given other variables (z's)
  reg_weights <- 1 - par_m
  Prec_upd <- crossprod(Bgam, sweep(Bgam, MARGIN = 1, STATS = reg_weights, FUN = "*"))
  E_reciproc_sigmasq <- par_sigmasq$shape / par_sigmasq$rate
  E_reciproc_tausq <- par_tausq$shape / par_tausq$rate
  E_lambda <-  (1/E_reciproc_sigmasq) * E_reciproc_tausq
  Prec_upd <- Prec_upd + E_lambda * Sgam
  mu_upd <- solve(Prec_upd, crossprod(Bgam, reg_weights * y))

  par_beta$Prec <- Prec_upd
  par_beta$mu <- mu_upd

  par_sigmasq$shape <- sum(reg_weights)/2
  par_sigmasq$rate <- sum(reg_weights * ( (y - Bgam %*% mu_upd)^2 + (1/E_reciproc_sigmasq) * colSums(t(Bgam) *solve(par_beta$Prec, t(Bgam)))))/2

  # [4] update tau^2
  par_tausq$rate <- (sum(Sgam * solve(Prec_upd)) * (1/E_reciproc_sigmasq) + sum(mu_upd * (Sgam %*% mu_upd)))/2
 }

 par_beta <- list(mu = par_beta$mu,  Sigma = solve(par_beta$Prec) * (1/E_reciproc_sigmasq))
 par_sigmasq <- list(shape = par_sigmasq$shape, rate = par_sigmasq$rate)
 par_tausq <- list(shape = par_tausq$shape, rate = par_tausq$rate)
 par_alpha <- list(a = par_pi$a, b = par_pi$b)

 # Calculate R^2
 yhat <- as.vector(Bgam %*% par_beta$mu)
 sst <- sum((y - mean(y))^2)
 ssr <- sum((y - yhat)^2)
 R2 <- 1 - (ssr / sst)

 VI_fit <- list(B = Bgam, par_coef = par_beta, par_sigmasq = par_sigmasq,
                par_tausq = par_tausq, par_alpha = par_alpha,
                match.prob = c(par_m), family = list(family = gam_naive$family$family,
                                             link = gam_naive$family$link),
                R2 = R2, call = list(call = match.call(),
                            formula = gamobj$formula))

 class(VI_fit) <- "spline_mixture"
 VI_fit
}
