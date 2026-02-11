#' GLM with Mixture-Based Linkage Error Adjustment
#'
#' Fits a generalized linear model (GLM) accounting for mismatch errors using a
#' mixture model framework in the secondary analysis setting. The variance-covariance
#' matrix is estimated using the sandwich formula.
#'
#' @param x Design matrix for the primary outcome model (numeric matrix or data frame).
#' @param y Response vector for the primary outcome model.
#' @param family A family object (e.g., \code{gaussian}, \code{binomial}) specifying
#'   the error distribution and link function. Can be a character string or a function.
#' @param z Design matrix for the mismatch indicator model (mismatch covariates).
#'   If NULL, an intercept-only model is assumed.
#' @param m.rate The assumed overall mismatch rate (a proportion between 0 and 1).
#'   If provided, it imposes a constraint on the mismatch model intercept.
#' @param safe.matches Logical vector; \code{TRUE} indicates a "safe match" (treated
#'   as definitely correct), \code{FALSE} indicates a potential mismatch.
#' @param control An optional list of control parameters. Arguments passed via \code{...}
#'   will override values in this list.
#'   \itemize{
#'     \item \code{max.iter}: Maximum EM iterations (default: 1000).
#'     \item \code{cmax.iter}: Maximum iterations for the subroutine in the constrained logistic regression function (default: 1000).
#'     \item \code{tol}: Convergence tolerance (default: 1e-4).
#'     \item \code{init.beta}: Initial parameter estimates for the outcome model.
#'     \item \code{init.gamma}: Initial parameter estimates for the mismatch indicator model.
#'     \item \code{fy}: Estimated marginal density of the response. If NULL, estimated using Kernel Density Estimation or parametric assumption.
#'   }
#' @param ... Additional arguments passed to \code{control}.
#'
#' @returns A list of results:
#' \item{coefficients}{A named vector of coefficients for the outcome model.}
#' \item{m.coefficients}{A named vector of coefficients for the mismatch indicator model (gamma).}
#' \item{match.prob}{The posterior correct match probabilities (weights) for each observation.}
#' \item{residuals}{The working residuals, defined as \code{y - fitted.values}.}
#' \item{fitted.values}{The fitted mean values of the outcome model, obtained by transforming the linear predictors by the inverse of the link function.}
#' \item{linear.predictors}{The linear fit on the link scale.}
#' \item{deviance}{The deviance of the weighted outcome model at convergence.}
#' \item{null.deviance}{The deviance of the weighted null outcome model.}
#' \item{var}{The estimated variance-covariance matrix of the parameters (sandwich estimator).}
#' \item{dispersion}{The estimated dispersion parameter (e.g., variance for Gaussian, 1/shape for Gamma).}
#' \item{objective}{A vector tracking the negative log pseudo-likelihood at each iteration of the EM algorithm.}
#' \item{converged}{Logical indicating if the EM algorithm converged within \code{max.iter}.}
#' \item{rank}{The numeric rank of the fitted linear model.}
#' \item{df.residual}{The residual degrees of freedom.}
#' \item{df.null}{The residual degrees of freedom for the null model.}
#' \item{family}{The \code{family} object used.}
#'
#' @references
#' Slawski, M.*, West, B. T., Bukke, P., Wang, Z., Diao, G., & 
#' Ben-David, E. (2025). A general framework for regression with mismatched
#' data based on mixture modelling. \emph{Journal of the Royal Statistical Society
#' Series A: Statistics in Society}, 188(3), 896-919. \doi{10.1093/jrsssa/qnae083}
#' 
#' Slawski, M.*, Diao, G., Ben-David, E. (2021). A pseudo-likelihood approach to
#' linear regression with partially shuffled data. \emph{Journal of Computational
#' and Graphical Statistics}. 30(4), 991-1003. \doi{10.1080/10618600.2020.1870482}
#'
#' @importFrom stats glm.fit dnorm dgamma dpois dbinom plogis qlogis density approxfun coef sd family
#' @importFrom stats gaussian binomial quasibinomial optim model.matrix model.frame model.response
#'
#' @examples
#' \dontrun{
#' # 1. Simulate Data
#' set.seed(123)
#' n <- 1000
#' x <- matrix(rnorm(n), ncol = 1)
#' true_beta <- c(1, 2) # Intercept=1, Slope=2
#' y <- 1 + 2 * x[,1] + rnorm(n, sd = 0.5)
#'
#' # 2. Introduce Linkage Errors (Shuffle 10% of data)
#' m_rate <- 0.10
#' idx_mismatch <- sample(1:n, size = n * m_rate)
#' y[idx_mismatch] <- sample(y[idx_mismatch]) # Shuffle y for mismatches
#'
#' # 3. Create Linkage Covariates (z)
#' # Assume we have a "matching score" z related to correctness
#' z <- matrix(rnorm(n, mean = 1), ncol = 1)
#'
#' # 4. Fit the Mixture GLM
#' fit <- glmMixture(x = cbind(1, x), y = y,
#'                   family = "gaussian",
#'                   z = cbind(1, z),
#'                   m.rate = 0.10)
#' }
#'
#' @export
glmMixture <- function(x, y, family,
                       z, m.rate = NULL, safe.matches = NULL,
                       control = list(), ...) {
  
  # Argument Handling and Setup
  dots <- list(...)
  # Merge (...) into control, giving (...) precedence
  con <- list(init.beta = NULL, init.gamma = NULL, fy = NULL,
              max.iter = 1000, tol = 1E-4, cmax.iter = 1000)
  
  # Update control with provided list, then override with dots
  con[names(control)] <- control
  con[names(dots)] <- dots
  
  x <- as.matrix(x)
  z <- as.matrix(z)
  y <- as.numeric(y)
  n <- nrow(x)
  p <- ncol(x)
  
  # Calculate logit bound for m.rate constraint
  if(!is.null(m.rate)){
    logitbound <- -log((1 - m.rate) / m.rate)
  } else {
    logitbound <- NULL
  }
  
  # Family Setup
  family_arg <- if (!is.null(family)) family else stats::gaussian()
  
  if (is.character(family_arg)) {
    if (tolower(family_arg) == "gamma") family <- stats::Gamma(link = "log")
    else {
      family_fun <- try(get(family_arg, mode = "function", envir = parent.frame()), silent = TRUE)
      if (inherits(family_fun, "try-error")) stop(paste("Invalid family argument:", family_arg))
      family <- family_fun()
    }
  } else if (is.function(family_arg)) {
    family <- family_arg()
  } else {
    family <- family_arg
  }
  if (is.null(family$family)) stop("Invalid family object.")
  
  if (is.null(safe.matches)) safe.matches <- rep(FALSE, n)
  
  # Marginal Density (fy) Estimation
  fy <- con$fy
  if (is.null(fy)) {
    if (family$family == "binomial") {
      fy <- stats::dbinom(y, size = 1, prob = mean(y))
    } else if (family$family == "gaussian") {
      fy <- stats::dnorm(y, mean = mean(y), sd = stats::sd(y))
    } else {
      # Use Kernel Density Estimation for generic cases
      kde_obj <- stats::density(y)
      approx_fun <- stats::approxfun(x = kde_obj$x, y = kde_obj$y, method = "linear", rule = 2)
      fy <- approx_fun(y)
      fy[fy < 1e-300] <- 1e-300
    }
  }
  
  # Initialization
  beta_cur <- con$init.beta
  dispersion_cur <- NULL
  shape_cur <- NULL
  
  if (is.null(beta_cur)) {
    # Standard GLM fit as starting point
    init_fit <- stats::glm.fit(x, y, family = family)
    beta_cur <- stats::coef(init_fit)
    mu_init <- init_fit$fitted.values
  } else {
    beta_cur <- as.vector(beta_cur)
    mu_init <- family$linkinv(x %*% beta_cur)
  }
  
  # Initialize dispersion/shape
  if (family$family == "gaussian") {
    dispersion_cur <- sum((y - mu_init)^2) / (n - p)
  } else if (family$family == "Gamma") {
    dispersion_est <- sum((y - mu_init)^2 / mu_init^2) / (n - p)
    shape_cur <- 1 / dispersion_est
  } else {
    dispersion_cur <- 1
  }
  
  # Initialize Gamma (Mismatch parameters)
  gamma_cur <- con$init.gamma
  if (is.null(gamma_cur)) {
    if (!is.null(logitbound)) {
      # If constraint exists, initialize close to bound
      if (ncol(z) == 1 && all(z == 1)) gamma_cur <- rep(-logitbound, ncol(z))
      else gamma_cur <- c(min(logitbound, 0), rep(0, ncol(z) - 1))
    } else {
      gamma_cur <- rep(0, ncol(z))
    }
  }
  
  # Helper functions
  h_func <- function(eta) {
    pi <- stats::plogis(eta)
    list(fun = pi, dfun = pi * (1 - pi), d2fun = pi * (1 - pi) * (1 - 2 * pi))
  }
  
  eval_phi <- function(eta, sub_idx, family_obj, sigma_sq = NULL, shape = NULL) {
    y_sub <- y[sub_idx]
    eta_sub <- eta[sub_idx]
    mu_sub <- family_obj$linkinv(eta_sub)
    
    if (family_obj$family == "gaussian") {
      fun <- stats::dnorm(y_sub, mean = mu_sub, sd = sqrt(sigma_sq))
    } else if (family_obj$family == "Gamma") {
      fun <- stats::dgamma(y_sub, shape = shape, rate = shape / mu_sub)
    } else if (family_obj$family == "poisson") {
      fun <- stats::dpois(y_sub, lambda = mu_sub)
    } else {
      fun <- stats::dbinom(y_sub, size = 1, prob = mu_sub)
    }
    
    mu_prime <- family_obj$mu.eta(eta_sub)
    var_mu <- family_obj$variance(mu_sub)
    
    if (family_obj$family == "Gamma") score_eta <- (y_sub - mu_sub) / var_mu * mu_prime * shape
    else score_eta <- (y_sub - mu_sub) / var_mu * mu_prime
    
    dfun <- fun * score_eta
    
    if (family_obj$family == "Gamma") expected_hess <- - (mu_prime^2) / var_mu * shape
    else expected_hess <- - (mu_prime^2) / var_mu
    
    d2fun <- fun * (score_eta^2 + expected_hess)
    res <- list(fun = fun, dfun = dfun, d2fun = d2fun)
    
    if (family_obj$family == "gaussian") {
      res$dfun_beta <- dfun / sigma_sq
      res$d2fun_beta <- d2fun / sigma_sq
      d_fun_sigma <- fun * ((y_sub - mu_sub)^2 / sigma_sq^1.5 - 1 / sqrt(sigma_sq))
      res$dfun_sigma <- d_fun_sigma
      res$d2fun_sigma <- d_fun_sigma * ((y_sub - mu_sub)^2 / sigma_sq^1.5 - 1 / sqrt(sigma_sq)) +
        fun * (1 / sigma_sq - 3 * (y_sub - mu_sub)^2 / sigma_sq^2)
      res$d2fun_beta_sigma <- fun * (-2 * (y_sub - mu_sub) / sigma_sq^1.5) +
        d_fun_sigma * (y_sub - mu_sub) / sigma_sq
    } else if (family_obj$family == "Gamma") {
      score_shape <- log(shape / mu_sub) + 1 - digamma(shape) + log(y_sub) - y_sub / mu_sub
      res$dfun_shape <- fun * score_shape
      res$d2fun_shape <- res$dfun_shape * score_shape + fun * (1/shape - trigamma(shape))
    }
    return(res)
  }
  
  calc_loglik <- function(mu, hs, sigma_sq = NULL, shape = NULL) {
    if (family$family == "gaussian") phi <- stats::dnorm(y, mean = mu, sd = sqrt(sigma_sq))
    else if (family$family == "poisson") phi <- stats::dpois(y, lambda = mu)
    else if (family$family == "binomial") phi <- stats::dbinom(y, size = 1, prob = mu)
    else if (family$family == "Gamma") phi <- stats::dgamma(y, shape = shape, rate = shape / mu)
    
    term_match <- hs[!safe.matches] * phi[!safe.matches]
    term_mismatch <- (1 - hs[!safe.matches]) * fy[!safe.matches]
    
    # Avoid log(0)
    lik_terms <- term_match + term_mismatch
    lik_terms[lik_terms <= 0] <- 1e-300
    
    sum(-log(lik_terms)) - sum(log(phi[safe.matches]))
  }
  
  # EM Algorithm
  eta_cur <- x %*% beta_cur
  mu_cur <- family$linkinv(eta_cur)
  
  # hs is P(m=0|z)
  hs <- h_func(z %*% gamma_cur)$fun
  hs[safe.matches] <- 1
  
  iter <- 1
  objs <- numeric(con$max.iter)
  p_cur <- rep(0, n) # Posterior probability of match
  w_glm_fit <- NULL
  
  if (family$family == "gaussian") objs[1] <- calc_loglik(mu_cur, hs, sigma_sq = dispersion_cur)
  else if (family$family == "Gamma") objs[1] <- calc_loglik(mu_cur, hs, shape = shape_cur)
  else objs[1] <- calc_loglik(mu_cur, hs)
  
  while (iter < con$max.iter) {
    # E-Step: Calculate posterior match probabilities (p_cur)
    if (family$family == "gaussian") phi <- stats::dnorm(y, mean = mu_cur, sd = sqrt(dispersion_cur))
    else if (family$family == "Gamma") phi <- stats::dgamma(y, shape = shape_cur, rate = shape_cur / mu_cur)
    else if (family$family == "poisson") phi <- stats::dpois(y, lambda = mu_cur)
    else phi <- stats::dbinom(y, size = 1, prob = mu_cur)
    
    num <- hs[!safe.matches] * phi[!safe.matches]
    denom <- num + (1 - hs[!safe.matches]) * fy[!safe.matches]
    ratio <- num / denom
    ratio[denom == 0] <- 0 # Handle division by zero
    
    p_cur[!safe.matches] <- ratio
    p_cur[safe.matches] <- 1
    
    if (anyNA(p_cur)) { warning("EM: NA weights."); break }
    
    # M-Step: Update Gamma (Mismatch Model)
    z_sub <- z[!safe.matches, , drop = FALSE]
    p_sub <- p_cur[!safe.matches]
    
    if (!is.null(logitbound)) {
      glm_h <- constrained_logistic_regression(z_sub, 1 - p_sub, logitbound, con$cmax.iter)
      # constrained_logistic_regression models mismatch (1-p), so coefs are for mismatch.
      # h_func models match. relation: gamma_match = -gamma_mismatch
      gamma_cur <- -glm_h$beta
    } else {
      # Use quasibinomial to handle fractional weights [0,1]
      glm_h <- stats::glm.fit(z_sub, p_sub, family = stats::quasibinomial())
      gamma_cur <- stats::coef(glm_h)
    }
    
    hs[!safe.matches] <- h_func(z[!safe.matches, , drop = FALSE] %*% gamma_cur)$fun
    
    # M-Step: Update Beta (Outcome Model)
    fit_fam <- if (family$family == "binomial") stats::quasibinomial(link = family$link) else family
    
    if (family$family == "Gamma") {
      w_glm_fit <- stats::glm.fit(x, y, weights = p_cur, family = fit_fam)
      beta_cur <- stats::coef(w_glm_fit)
      eta_cur <- x %*% beta_cur
      mu_cur <- family$linkinv(eta_cur)
      
      # Update Shape for Gamma
      sum_w <- sum(p_cur)
      term_k <- sum(p_cur * (log(y / mu_cur) - y / mu_cur))
      
      if (is.null(shape_cur) || is.na(shape_cur) || shape_cur <= 0) {
        shape_cur <- 1 / (sum((y - mu_cur)^2 / mu_cur^2) / (n - p))
      }
      
      # Newton-Raphson for Shape
      for (i in 1:20) {
        sc <- sum_w * (log(shape_cur) + 1 - digamma(shape_cur)) + term_k
        hess <- sum_w * (1 / shape_cur - trigamma(shape_cur))
        new_s <- shape_cur - sc / hess
        sf <- 1
        while (new_s <= 1e-4 && sf > 1e-4) { sf <- sf * 0.5; new_s <- shape_cur - sc / hess * sf }
        if (new_s <= 1e-4) new_s <- 1e-4
        if (abs(new_s - shape_cur) < 1e-6) { shape_cur <- new_s; break }
        shape_cur <- new_s
      }
    } else {
      w_glm_fit <- stats::glm.fit(x, y, weights = p_cur, family = fit_fam)
      beta_cur <- stats::coef(w_glm_fit)
      eta_cur <- x %*% beta_cur
      mu_cur <- family$linkinv(eta_cur)
    }
    
    if (family$family == "gaussian") {
      # Update Dispersion (Variance)
      dispersion_cur <- sum(p_cur * (y - mu_cur)^2) / sum(p_cur)
    }
    
    iter <- iter + 1
    
    # Calculate Objective
    if (family$family == "gaussian") objs[iter] <- calc_loglik(mu_cur, hs, sigma_sq = dispersion_cur)
    else if (family$family == "Gamma") objs[iter] <- calc_loglik(mu_cur, hs, shape = shape_cur)
    else objs[iter] <- calc_loglik(mu_cur, hs)
    
    if (!is.na(objs[iter]) && abs(objs[iter] - objs[iter - 1]) < con$tol) break
  }
  
  # Variance-Covariance Estimation (Sandwich)
  evals <- eval_phi(eta_cur, !safe.matches, family,
                    sigma_sq = if (family$family == "gaussian") dispersion_cur else NULL,
                    shape = shape_cur)
  
  # Note: drop=FALSE is required in all subsetting below to handle single-mismatch cases
  z_sub <- z[!safe.matches, , drop = FALSE] 
  x_sub <- x[!safe.matches, , drop = FALSE]
  
  h_eval <- h_func(z_sub %*% gamma_cur)
  mix_prob <- fy[!safe.matches] * (1 - h_eval$fun) + h_eval$fun * evals$fun
  
  if (family$family == "gaussian") {
    w_beta <- (-1) * evals$dfun_beta * h_eval$fun / mix_prob
    w_extra <- (-1) * evals$dfun_sigma * h_eval$fun / mix_prob
  } else if (family$family == "Gamma") {
    w_beta <- (-1) * evals$dfun * h_eval$fun / mix_prob
    w_extra <- (-1) * evals$dfun_shape * h_eval$fun / mix_prob
  } else {
    w_beta <- (-1) * evals$dfun * h_eval$fun / mix_prob
  }
  
  # Gradient w.r.t gamma
  w_gamma <- ((-1) * (evals$fun - fy[!safe.matches]) * h_eval$dfun) / mix_prob
  
  x_w1 <- sweep(x_sub, 1, w_beta, "*")
  Delta_w3 <- sweep(z_sub, 1, w_gamma, "*")
  
  if (family$family %in% c("gaussian", "Gamma")) {
    x_w2 <- matrix(w_extra, ncol = 1)
    meat <- crossprod(cbind(x_w1, x_w2, Delta_w3))
  } else {
    meat <- crossprod(cbind(x_w1, Delta_w3))
  }
  
  # Second derivatives for Hessian
  if (family$family == "gaussian") {
    w_beta2 <- (-(h_eval$fun * evals$d2fun_beta) / mix_prob) + w_beta^2
    w_extra2 <- (-(h_eval$fun * evals$d2fun_sigma) / mix_prob) + w_extra^2
    w_beta_extra <- (-(evals$d2fun_beta_sigma * h_eval$fun) / mix_prob) + w_beta * w_extra
  } else if (family$family == "Gamma") {
    w_beta2 <- (-(h_eval$fun * evals$d2fun) / mix_prob) + w_beta^2
    w_extra2 <- (-(h_eval$fun * evals$d2fun_shape) / mix_prob) + w_extra^2
    w_beta_extra <- w_beta * w_extra
  } else {
    w_beta2 <- (-(h_eval$fun * evals$d2fun) / mix_prob) + w_beta^2
  }
  w_gamma2 <- (-(evals$fun - fy[!safe.matches]) * h_eval$d2fun / mix_prob) + w_gamma^2
  
  if (family$family == "gaussian") term_bg <- evals$dfun_beta
  else term_bg <- evals$dfun
  
  w_beta_gamma <- (-(term_bg * h_eval$dfun) / mix_prob) +
    ((evals$fun - fy[!safe.matches]) * h_eval$fun * h_eval$dfun * term_bg) / (mix_prob^2)
  
  x_w4 <- sweep(x_sub, 1, w_beta2, "*")
  Delta_w6 <- sweep(z_sub, 1, w_gamma2, "*")
  x_w5 <- sweep(x_sub, 1, w_beta_gamma, "*")
  
  d_x <- ncol(x); d_z <- ncol(z)
  
  # Construct Hessian (Use x_sub and z_sub instead of subsetting inline)
  if (family$family %in% c("gaussian", "Gamma")) {
    size <- d_x + 1 + d_z
    Hess <- matrix(0, size, size)
    one <- matrix(1, sum(!safe.matches), 1)
    
    Hess[1:d_x, 1:d_x] <- crossprod(x_sub, x_w4)
    Hess[d_x + 1, d_x + 1] <- crossprod(one, one * w_extra2)
    Hess[(d_x + 2):size, (d_x + 2):size] <- crossprod(z_sub, Delta_w6)
    
    bg <- crossprod(x_w5, z_sub)
    Hess[1:d_x, (d_x + 2):size] <- bg
    Hess[(d_x + 2):size, 1:d_x] <- t(bg)
    
    bs <- crossprod(x_sub, one * w_beta_extra)
    Hess[1:d_x, d_x + 1] <- bs
    Hess[d_x + 1, 1:d_x] <- t(bs)
    
    sg <- crossprod(one * (w_extra * w_gamma), z_sub)
    Hess[d_x + 1, (d_x + 2):size] <- sg
    Hess[(d_x + 2):size, d_x + 1] <- t(sg)
    
  } else {
    size <- d_x + d_z
    Hess <- matrix(0, size, size)
    Hess[1:d_x, 1:d_x] <- crossprod(x_sub, x_w4)
    Hess[(d_x + 1):size, (d_x + 1):size] <- crossprod(z_sub, Delta_w6)
    
    bg <- crossprod(x_w5, z_sub)
    Hess[1:d_x, (d_x + 1):size] <- bg
    Hess[(d_x + 1):size, 1:d_x] <- t(bg)
  }
  
  cov_1_hat <- try(solve(Hess, meat), silent = TRUE)
  covhat <- if (inherits(cov_1_hat, "try-error")) matrix(NA, nrow(Hess), ncol(Hess)) else t(solve(Hess, t(cov_1_hat)))
  
  beta_n <- colnames(x); if (is.null(beta_n)) beta_n <- paste0("beta", 1:p)
  gamma_n <- colnames(z); if (is.null(gamma_n)) gamma_n <- paste0("gamma", 1:ncol(z))
  
  if (family$family == "gaussian") {
    idx <- p + 1; sc <- 2 * sqrt(dispersion_cur)
    if (!anyNA(covhat)) { covhat[idx, ] <- covhat[idx, ] * sc; covhat[, idx] <- covhat[, idx] * sc }
    rn <- c(paste("coef", beta_n), "dispersion", paste("m.coef", gamma_n))
  } else if (family$family == "Gamma") {
    rn <- c(paste("coef", beta_n), "shape", paste("m.coef", gamma_n))
  } else {
    rn <- c(paste("coef", beta_n), paste("m.coef", gamma_n))
  }
  rownames(covhat) <- colnames(covhat) <- rn
  
  # Deviance of the outcome model: Use the weighted deviance from the M-step.
  #    This reflects the fit of the outcome model conditional on the posterior match probabilities.
  deviance_model <- w_glm_fit$deviance
  df_residual <- w_glm_fit$df.residual
  
  # Null Deviance: Fit weighted null model
  #    Standard glm null deviance is based on intercept-only model with same weights
  null_deviance <- w_glm_fit$null.deviance
  df_null <- w_glm_fit$df.null
  
  out <- list(coefficients = beta_cur, 
              m.coefficients = gamma_cur, 
              residuals = y - mu_cur,
              linear.predictors = eta_cur,
              fitted.values = mu_cur, 
              deviance = deviance_model,
              null.deviance = null_deviance,            
              df.residual = n - p, 
              df.null = n - 1,
              rank = p, 
              family = family, 
              converged = iter < con$max.iter, 
              match.prob = hs, 
              var = covhat,
              objective = objs[1:iter],
              call = match.call())
  
  if (family$family == "gaussian") out$dispersion <- dispersion_cur
  else if (family$family == "Gamma") out$dispersion <- 1 / shape_cur
  else out$dispersion <- 1
  
  # Assign class inheritance
  class(out) <- c("glmMixture")
  
  return(out)
}

#' @keywords internal
#' @export
fitglm.adjMixture <- function(x, y, family, adjustment, control, ...) {
  
  # -------------------------------------------------------------------------
  # 1. Validation and Data Retrieval
  # -------------------------------------------------------------------------
  full_data <- adjustment$data_ref$data
  if (is.null(full_data)) {
    stop("The 'adjustment' object does not contain linked data. ",
         "Please recreate the object with 'linked.data' provided.", call. = FALSE)
  }
  
  # -------------------------------------------------------------------------
  # 2. Stage 1: Align Adjustment Data to Outcome Model (X, Y)
  # -------------------------------------------------------------------------
  # plglm() has already applied 'subset' and 'na.action' to x and y.
  # We use row names to synchronize the adjustment data.
  
  subset_names <- rownames(x)
  
  # If x has no row names (rare), assume 1:1 mapping if sizes match.
  # This handles cases where users strip names or use matrices without names.
  if (is.null(subset_names)) {
    if (nrow(x) != nrow(full_data)) {
      stop("Row mismatch: Model matrix 'x' has no row names and its length (", nrow(x), 
           ") differs from the adjustment data (", nrow(full_data), "). ",
           "Ensure 'linked.data' matches the data passed to 'plglm'.", call. = FALSE)
    }
    # Assume the user provided the exact same dataset in the same order
    data_subset <- full_data
  } else {
    # Match by name. Using strict matching ensures order is preserved.
    # Note: data.frames are guaranteed to have unique row names in R.
    idx_map <- match(subset_names, rownames(full_data))
    
    if (anyNA(idx_map)) {
      stop("Row mismatch: Some observations in the model matrix could not be matched ",
           "to the adjustment data. This usually happens if 'data' in plglm() ",
           "is different from the data used to create the adjustment object.", call. = FALSE)
    }
    data_subset <- full_data[idx_map, , drop = FALSE]
  }
  
  # -------------------------------------------------------------------------
  # 3. Construct Mismatch Covariates (Z)
  # -------------------------------------------------------------------------
  # We use model.frame/model.matrix to handle factors and transformations in m.formula
  m_formula <- adjustment$m.formula
  
  # We use 'na.pass' here because we want to detect NAs manually in the next step
  # to ensure we drop the corresponding rows in X and Y simultaneously.
  Z_frame <- tryCatch({
    stats::model.frame(m_formula, data = data_subset, na.action = stats::na.pass)
  }, error = function(e) {
    stop("Failed to resolve variables in 'm.formula': ", e$message, call. = FALSE)
  })
  
  Z <- stats::model.matrix(m_formula, Z_frame)
  
  # -------------------------------------------------------------------------
  # 4. Stage 2: Secondary Intersection (Handle Missingness in Z)
  # -------------------------------------------------------------------------
  # X and Y are already clean (no NAs). But Z might have NAs.
  # If Z has NAs, we must drop those rows from X, Y, and Z to stay aligned.
  
  # Identify complete cases in Z
  keep_idx <- stats::complete.cases(Z)
  
  if (!all(keep_idx)) {
    n_dropped <- sum(!keep_idx)
    warning(sprintf("Dropped %d observation(s) due to missing values in the 
                    mismatch covariates (Z) that were not missing in the 
                    outcome variables.", n_dropped), call. = FALSE)
    
    # Apply the subset
    x <- x[keep_idx, , drop = FALSE]
    y <- y[keep_idx]
    Z <- Z[keep_idx, , drop = FALSE]
    
    # Also update the row map for the safe matches step below
    if (!is.null(subset_names)) {
      subset_names <- subset_names[keep_idx]
      idx_map <- idx_map[keep_idx]
    } else {
      # If we were using implicit ordering, we must subset the implicit map
      # This logic handles the "no row names" edge case
      full_indices <- seq_len(nrow(full_data))
      idx_map <- full_indices[keep_idx] 
    }
  }
  
  # -------------------------------------------------------------------------
  # 5. Process Safe Matches
  # -------------------------------------------------------------------------
  safe_matches_all <- adjustment$safe.matches
  safe_matches_sub <- NULL
  
  if (!is.null(safe_matches_all)) {
    # Ensure the safe.matches vector aligns with the final subset of data
    if (!is.null(subset_names)) {
      # We use the updated idx_map which accounts for both plglm subsetting AND Z-missingness
      safe_matches_sub <- safe_matches_all[idx_map]
    } else {
      # Fallback for no row names
      safe_matches_sub <- safe_matches_all[1:nrow(x)]
    }
  }
  
  # -------------------------------------------------------------------------
  # 6. Dispatch to Computational Engine
  # -------------------------------------------------------------------------
  fit <- glmMixture(
    x = x,
    y = y,
    family = family,
    z = Z,
    m.rate = adjustment$m.rate,
    safe.matches = safe_matches_sub,
    control = control,
    ...
  )
  
  # -------------------------------------------------------------------------
  # 7. Post-Processing
  # -------------------------------------------------------------------------
  fit$adjustment <- adjustment
  fit$m.formula  <- m_formula
  fit$call <- match.call() 
  
  # Restore row names to observation-level outputs
  if (!is.null(subset_names)) {
    names(fit$residuals)         <- subset_names
    names(fit$fitted.values)     <- subset_names
    names(fit$linear.predictors) <- subset_names
    # match.prob is a vector in glmMixture, verify before naming
    if (!is.null(fit$match.prob)) names(fit$match.prob) <- subset_names
  }
  
  return(fit)
}