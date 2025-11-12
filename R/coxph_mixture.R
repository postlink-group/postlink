#' Fit a CoxPH Model Using a Mixture-Modeling Approach
#' @import survival
#' @description
#' Fit Cox proportional hazards regression adjusted for mismatched data.
#' Information about the underlying record linkage process can be
#' incorporated into the method if available (e.g., assumed overall mismatch rate,
#' safe matches, predictors of match status, or predicted probabilities of correct
#' matches).
#'
#' @param formula a formula object for the outcome model, the response should be
#' provided using the `Surv` function and the covariates should be separated by + signs.
#' @param data a data.frame with linked data used in "formula" and "formula.m" (optional)
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
#' @returns a list of results from the function called depending on the "family"
#' specified.
#' \item{coefficients}{the outcome model coefficient estimates}
#' \item{match.prob}{the posterior correct match probabilities for observations
#' given parameter estimates}
#' \item{objective}{a variable that tracks the negative log pseudo-likelihood
#' for all iterations of the EM algorithm.}
#' \item{standard.errors}{the estimated standard errors}
#' \item{m.coefficients}{the correct match model coefficient estimates}
#' \item{call}{the matched call}
#' \item{wfit}{an internal-use object for the predict function}
#' \item{vcov}{the variance-covariance matrix}
#' \item{Lambdahat_0}{the baseline cumulative hazard (using weighted Breslow estimator) when the family is "cox"}
#' \item{g_Lambdahat_0}{the baseline cumulative hazard for the marginal density
#' of the response variable (using Nelson-Aalen estimator)}
#'
#' @note
#' The reference below discuss the implemented framework in more detail. The standard
#' errors are estimated Louis' method for the "cox" family (Bukke et al., 2023).\cr
#'
#' @references Bukke, P., Ben-David, E., Diao, G., Slawski, M., & West, B. T. (2025).
#' Cox Proportional Hazards Regression Using Linked Data: An Approach Based on Mixture Modelling.\cr
#'
#' @export
coxph_mixture <- function(formula, data,
                          m.formula, safe.matches, m.rate,
                          control = list(init.beta = "default",
                                         init.gamma = "default",
                                         fy = "default",
                                         max.iter = 1000,
                                         tol = 1E-4, cmax.iter = 1000),...){
 # 1. SET-UP
 # -----------------------------------------------------------------------------
 # CONTROLS
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

 # MAIN ARGUMENTS
 if(!missing(m.rate) && (m.rate <= 0 | m.rate >= 1)){
  stop("Error: assumed mismatch rate should be a proportion between 0 and 1")}

 if(missing(formula)){ stop("Error: a formula for the outcome model is required")}
 if(!inherits(formula, "formula")){ stop("Error: formula should be a formula object")}

 if(!missing(data) && (!is.data.frame(data) & !is.list(data))){
  stop("Error: data should be a data.frame or list")}

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

 # INPUTS
 # formula, data, and m.formula (X, Y, cens, logis_ps, n, p)
  mf <- model.frame(formula, data = data)
  cens <- 1 - mf[[1]][,"status"]
  X <- as.matrix(mf[,-1])
  y <- mf[[1]][,"time"]
  if(attr(mf[[1]], "type") != "right"){
   stop("Error: censoring type other than right-censoring is
        not currently supported")
  }
  if(!is.Surv(mf[[1]])){
   stop(("Error: response should be a survival object"))}
  n <- nrow(X)
  p <- ncol(X)

  if(missing(m.formula)){
   logis_ps <- matrix(nrow = n, ncol = 1, data = 1)
   colnames(logis_ps) <- "(Intercept)"
  } else {
   if(!is.null(model.response(model.frame(m.formula, data = data)))){
    stop("Error: m.formula should be a one-sided formula")}
   logis_ps <- model.matrix(m.formula, data=data)}

 if(any(is.na(X)) | any(is.na(y))){stop("Error (formula): Cannot have
  missing observations")}
 if(any(is.na(logis_ps))){stop("Error (m.formula): Cannot have missing observations")}
 if(nrow(logis_ps) != n){stop("Error (m.formula): Number of observations in formula
  and m.formula data should match")}

 # safe matches (is_flagged)
  if(missing(safe.matches)){
   is_flagged <- rep(FALSE,n)
 } else {
    is_flagged <- safe.matches
  }

 if(any(is.na(is_flagged))){stop("Error (safe.matches): Cannot have
  missing observations")}
 if(length(is_flagged) != n){stop("Error (safe.matches): Length of safe.matches
  should match number of observations")}

 # m.rate (logitbound)
 if(!missing(m.rate)){
  logitbound <- -log((1 - m.rate)/m.rate)
 }

 # fy
 if (!identical(fy, "default")){
  fy <- as.vector(fy)
  if(length(fy) != n){
   warning("Default 'fy' used. 'fy' should be a vector with a value
              for each observation")
   fy <- "default"}
 }
 if (identical(fy, "default")){
  survf0 <- survfit(Surv(y, event = 1 - cens) ~ 1)
  times <- survf0$time
  matchytimes <- sapply(y, function(x) which.min(abs(times - x)))
  D <- c(survf0$cumhaz[1], diff(survf0$cumhaz)/diff(times))

  g_lambdahat_0 <- D[matchytimes]
  g_Lambdahat_0 <- survf0$cumhaz[matchytimes]

  fy <- g_lambdahat_0^(1-cens) * exp(-g_Lambdahat_0)
 }

 # init.beta
 if (!identical(init.beta, "default")){
  betacur <- as.vector(init.beta)
  if(length(betacur) != p){
   warning("Default 'init.beta' used. 'init.beta' should be a vector
           of length ", p)
   init.beta <- "default"}
 }
 if (identical(init.beta, "default")){
  creg <- coxph(Surv(y, event = 1 - cens) ~ X)
  beta_cur <- creg$coef
 }

 # init.gamma
 Delta <- logis_ps
 if (!identical(init.gamma, "default")){
  gammacur <- as.vector(init.gamma)
  if(length(gammacur) != p){
   warning("Default 'init.gamma' used. 'init.gamma' should be a vector
           of length ", ncol(Delta))
   init.gamma <- "default"}
 }
 if (identical(init.gamma, "default")){
  if(!missing(m.rate)){
   if(identical(Delta, matrix(nrow = n, ncol = 1, data = 1))){
    gammacur <- rep(-logitbound, ncol(Delta))
   } else {
    gammacur <- c(min(logitbound, 0), rep(0, ncol(Delta)-1))
   }
  } else {
   gammacur <- rep(0, ncol(Delta))
  }
 }

 # 2. INITIALIZE EM-ALGORITHM
 # -------------------------------------------------------------------------
 fymu <- function(mu, cens, l0, L0) (exp(mu) * l0)^(1 - cens) *
  exp(-exp(mu) * L0)

 nloglik <- function(mu, cens, hs, l0, L0) sum(-log(hs[!is_flagged] *
                                                     fymu(mu, cens, l0, L0)[!is_flagged] + (1 - hs[!is_flagged])*
                                                     fy[!is_flagged])) - sum(log(fymu(mu, cens, l0, L0)[is_flagged]))

 ## hs
 hgamma <- function(eta){
  pi <- plogis(eta)
  return(list(fun = pi, dfun = pi*(1-pi), d2fun = pi*(1-pi)*(1 - 2*pi)))
 }
 hs <- hgamma(Delta %*% gammacur)$fun
 hs[is_flagged] <- 1

 ## mucur
 mu <- X %*% beta_cur

 ## objs
 iter <- 1
 Breslow_Estimator <- basehaz(creg, centered=FALSE)
 cumhazard <- Breslow_Estimator$hazard
 times <- Breslow_Estimator$time
 matchytimes <- sapply(y, function(x) which.min(abs(times - x)))
 D1 <-c(cumhazard[1], diff(cumhazard)/diff(times))
 lambdahat_0 <- D1[matchytimes]
 Lambdahat_0 <- cumhazard[matchytimes]
 nloglik_cur <- nloglik(mu, cens, hs, lambdahat_0, Lambdahat_0)
 objs <- numeric(max.iter)
 objs[iter] <- nloglik_cur

 ## pcur
 pcur = rep(0,n)

 ## additional trackers
 track_beta <- matrix(0, nrow = max.iter, ncol = p)
 track_lam <- matrix(0, nrow = max.iter, ncol = n)
 track_Lam <- matrix(0, nrow = max.iter, ncol = n)

 track_beta[1,] <- beta_cur
 track_lam[1,] <- lambdahat_0
 track_Lam[1,] <- Lambdahat_0

 lambdahat_0_ <- lambdahat_0
 Lambdahat_0_ <- Lambdahat_0

 # 3. EM-ALGORITHM
 # -------------------------------------------------------------------------
 while(iter < max.iter){
  num <- hs[!is_flagged] * fymu(mu, cens, lambdahat_0_,
                                Lambdahat_0_)[!is_flagged]
  denom <- num + (1-hs[!is_flagged]) * fy[!is_flagged]
  pcur[!is_flagged] <- num/denom
  pcur[is_flagged] <- 1

  if(anyNA(pcur)){
   warning("EM algorithm did not converge. NA observation weight(s) occurred
           in the E-step.")
   return(list(coefficients =  beta_cur, m.coefficients = gammacur,
               match.prob = hs, objective = objs[1:(iter)],
               Lambdahat_0 = Lambdahat_0_,  g_Lambdahat_0= g_Lambdahat_0))
  }

  if(!missing(m.rate)){
   if(identical(Delta, matrix(nrow = n, ncol = 1, data = 1))){
    glm_h <- glm(pcur[!is_flagged] ~ Delta[!is_flagged,] - 1,
                 family = quasibinomial)
    gammacur <- max(coef(glm_h), -logitbound)
   }
   else {
    glm_h <- constrained_logistic_regression(Delta[!is_flagged,],
                                             1-pcur[!is_flagged], logitbound,
                                             cmax.iter)
    gammacur <- -glm_h$beta
   }
  } else {
   glm_h <- glm(pcur[!is_flagged] ~ Delta[!is_flagged,] - 1,
                family = quasibinomial)
   gammacur <- coef(glm_h)
  }

  hs[!is_flagged] <- hgamma(Delta[!is_flagged,,drop=FALSE] %*% gammacur)$fun

  creg <- coxph(Surv(y, event = 1 - cens) ~ X, weights =
                 pmax(drop(pcur),1E-6), x = TRUE)
  mu <- X %*% creg$coefficients
  beta_cur <- creg$coefficients
  track_beta[iter+1,] <- beta_cur

  Breslow_Estimator <- basehaz(creg, centered=FALSE)
  cumhazard <- Breslow_Estimator$hazard
  times <- Breslow_Estimator$time
  matchytimes <- sapply(y, function(x) which.min(abs(times - x)))
  D1 <- c(cumhazard[1], diff(cumhazard)/diff(times))

  lambdahat_0_ <- D1[matchytimes]
  Lambdahat_0_ <- cumhazard[matchytimes]

  track_lam[iter+1,] <- lambdahat_0_
  track_Lam[iter+1,] <- Lambdahat_0_

  iter <- iter + 1
  objs[iter] <- nloglik(mu, cens, hs, lambdahat_0_, Lambdahat_0_)

  if(is.na(objs[iter])){
   warning("EM algorithm did not converge. NA objective value occurred.")
   return(list(coefficients =  beta_cur, m.coefficients = gammacur,
               match.prob = hs, objective = objs[1:(iter)],
               Lambdahat_0 = Lambdahat_0_,  g_Lambdahat_0= g_Lambdahat_0))
  }

  if(objs[iter] + tol > objs[iter-1]){
   break
  }

 }
 names(gammacur) <- colnames(logis_ps)

 # 4. STANDARD ERRORS
 # -------------------------------------------------------------------------
 Xorig <- X # design matrix
 yorig <- y # response
 censoring <- cens # censoring indicator (1 if yes, 0 if no)
 muorig <- mu # Xorig %*% betacur
 Z <- Delta # m-model covariates

 ## PART 1: Run k Monte-Carlo Samples (i.e., of the {m_i}_{i=1}^n's)
 k <- 1000

 for (sample in 1:k){
  # 1. Sample {m_i^[k]}_{i=1}^n given data, \hat{\beta}, \hat{\gamma}
  # 1-pcur is conditional probability of mismatch
  m_i <- rbinom(n, size = 1, prob = (1-pcur))

  # 2. Subset Data w/ m_i = 0
  X <- Xorig[m_i == 0,,drop = FALSE]
  y <- yorig[m_i == 0,drop = FALSE]
  cens <- censoring[m_i == 0,drop = FALSE]
  mu <- muorig[m_i == 0,,drop = FALSE]

  zgam <- Z %*% gammacur

  pb <- length(beta_cur)
  pg  <- length(gammacur)
  pt <- pb + pg

  delta <- 1-cens # event indicator (1 if yes, 0 if no)
  n_event <- sum(delta)
  risk_set <- apply(as.matrix(y[delta == 1]), 1, function(yi)
   which(y >= yi)) # find R(y_i)

  # 3. Evaluate Gradient w.r.t \hat{\beta} & \hat{\gamma}
  #   derivative of -ve partial log-likelihood
  s_bg <- numeric(pt)
  val2 <- numeric(n_event)
  for (d1 in 1:pt){
   if (d1 <= pb){
    for(i in 1:n_event){
     val2[i] <- sum(exp(mu[risk_set[[i]]])*X[risk_set[[i]],d1])/
      sum(exp(mu[risk_set[[i]]]))
    }
    s_bg[d1] <- -sum(X[(delta ==1),d1]) + sum(val2)
   } else {
    indz <- d1 - pb
    s_bg[d1] <- sum((exp(zgam)*Z[,indz])/(1+exp(zgam)) - (1-m_i)*Z[,indz])
   }
  }
  if (sample == 1){
   gradient <- s_bg
  } else{
   gradient <- rbind(gradient, s_bg)
  }

  # 4. Evaluate Hessian w.r.t \hat{\beta} & \hat{\gamma}
  #    second derivative of -ve partial log-likelihood
  h_bg <- matrix(NA, nrow = pt, ncol = pt)
  low <- numeric(n_event)
  high <- numeric(n_event)
  dhigh <- numeric(n_event)
  dlow <- numeric(n_event)
  val <- numeric(n_event)
  for (d1 in 1:pt){
   for (d2 in 1:pt){
    if (d1 <= pb & d2 <= pb){
     for(i in 1:n_event){
      low[i] <- sum(exp(mu[risk_set[[i]]]))
      dlow[i] <- sum(exp(mu[risk_set[[i]]])*X[risk_set[[i]],d2])
      high[i] <- sum(exp(mu[risk_set[[i]]])*X[risk_set[[i]],d1])
      dhigh[i] <- sum(exp(mu[risk_set[[i]]])*X[risk_set[[i]],d1]*
                       X[risk_set[[i]],d2])
      val[i] <- (low[i]*dhigh[i] - high[i]*dlow[i])/(low[i])^2
     }
     h_bg[d1,d2] <- sum(val)
    }

    if (d1 <= pb & d2 > pb | d1 > pb & d2 <= pb){
     h_bg[d1,d2] <- 0
    }

    if (d1 > pb & d2 > pb){
     indz1 <- d1 - pb
     indz2 <- d2 - pb
     low <- (1+exp(zgam))
     dlow <- (exp(zgam)*Z[,indz2])
     high <- (exp(zgam)*Z[,indz1])
     dhigh <- (exp(zgam)*Z[,indz1]*Z[,indz2])
     h_bg[d1,d2] <- sum((low *dhigh  - high *dlow)/(low)^2)
    }
   }
  }

  if (sample == 1){
   hessian <- h_bg
  } else{
   hessian <- array(c(hessian, h_bg), dim = c(pt, pt, sample))
  }
 }
 row.names(gradient) <- NULL
 hmcg <- apply(hessian,c(1,2),mean) - cov(gradient)
 covhat <- solve(t(hmcg))
 se <- sqrt(diag(covhat))

 rownames(covhat) <- c(paste("coef", names(beta_cur)), paste("m.coef", names(gammacur)))
 colnames(covhat) <- c(paste("coef", names(beta_cur)), paste("m.coef", names(gammacur)))

 # 5. OUTPUTS
 # -------------------------------------------------------------------------
 x <- list(coefficients =  beta_cur, m.coefficients = gammacur,
           match.prob = c(hs), objective = objs[1:(iter)],
           Lambdahat_0 = Lambdahat_0_,  g_Lambdahat_0= g_Lambdahat_0,
           standard.errors = se, vcov = covhat,
           wfit = creg)
 x <- append(x, match.call())
 names(x)[[length(x)]] <- "call"

 class(x) <- "coxph_mixture"
 x

}
