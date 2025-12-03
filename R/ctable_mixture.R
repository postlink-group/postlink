#' Contingency Table Analysis Based on Mixture-Modeling
#' @description
#' Estimate the cell probabilities p_ij underlying a two-way contingency table under a multinomial sampling scheme
#' (i.e., the total number of counts is considered fixed). Note that in the absence of mismatches, this corresponds
#' to fitting a saturated model (the estimated cell probabilities will be exactly equal to the corresponding relative
#' frequencies). This estimator will be consistent under regularity conditions for correctly matched data. In the
#' presence of mismatches, however, the relative frequencies in the table will be biased. An adjustment is made by
#' fitting a mixture between a saturated model and an independence model. The mixture proportions for these models are
#' "1 – m.rate" and "m.rate", respectively. The mismatch rate cannot be estimated but needs to be pre-specified.
#'
#' @param formula a formula object with the left and right hand sides specifying the column
#' and row variable of the flat table, respectively.
#' @param data a data frame, list, environment, or contingency table (see `ftable.formula`).
#' @param m.rate the assumed mismatch rate. It should be greater than 0 and less than 1.
#' @param control an optional list variable for the EM algorithm to customize the maximum
#' iterations ("max.iter") and convergence tolerance for termination ("tol").
#' @param ... further arguments to the default ftable method may also be passed as arguments,
#' see `ftable.formula` and `ftable.default`. Also, control settings may be passed directly
#' instead of through the control argument.
#'
#' @returns a list of results.
#' \item{phat}{a matrix of the estimated cell probabilities.}
#' \item{phat0}{a matrix of the estimated cell probabilities for mismatched pairs under
#' an independence model. The corresponding parameter is a nuisance parameter.}
#' \item{vcov_phat}{the variance-covariance matrix of phat.}
#' \item{ftable}{a flat contingency table in matrix format containing the estimated counts
#' of each combination of the variables’ levels.}
#' \item{objs}{the sequence of objective values for the non-negative likelihood, which is
#' expected to be decreasing.}
#' \item{m.rate}{the mismatch rate assumed.}
#' \item{call}{the matched call.}
#'
#' @note
#' The reference below discusses the underlying methodology in more detail.
#'
#' @references Slawski, M., West, B. T., Bukke, P., Diao, G., Wang, Z., & Ben-David, E. (2023).
#' A General Framework for Regression with Mismatched Data Based on Mixture Modeling.
#' Under Review. < \doi{10.48550/arXiv.2306.00909} >\cr
#'
#' @export
ctable_mixture <- function(formula, data = NULL, m.rate,
                           control = list(max.iter = 100, tol = 1E-4),...){
 dcontrols <- list(...)
 if ("max.iter" %in% names(dcontrols)){
  max.iter <- dcontrols$max.iter
 } else {
  max.iter <- ifelse("max.iter" %in% names(control), control$max.iter, 1000)
 }

 if(max.iter < 2 | (floor(max.iter) != max.iter)){
  warning("Default max.iter used. max.iter should be an integer greater than 2")
  max.iter <- 1000}

 if ("tol" %in% names(dcontrols)){
  tol <- dcontrols$tol
 } else {
  tol <- ifelse("tol" %in% names(control), control$tol, 1E-4)
 }

 if(missing(m.rate)){
  stop("Error: the assumed mismatch rate is required")}

 if(!missing(m.rate) && (m.rate <= 0 | m.rate >= 1)){
  stop("Error: the assumed mismatch rate should be a proportion between 0 and 1")}

 if(missing(formula)){ stop("Error: a formula for the outcome model is required")}
 if(!inherits(formula, "formula")){ stop("Error: formula should be a formula object")}

 ctab <- if(!missing(data)){ftable(formula, data,...)} else{ftable(formula,...)}
 if(length(formula[[1]]) != 1 | length(formula[[2]]) != 1){
  stop("Error: the formula should contain one row variable and one column variable")
 }

 K <- nrow(ctab)
 L <- ncol(ctab)
 counts <- c(t(ctab))
 n <- sum(counts)

 alpha <- m.rate

 design <- expand.grid(1:L, 1:K)[,c(2,1)]
 design <- apply(design, 2, as.factor)
 colnames(design) <- c("x", "y")
 design <- as.data.frame(design)
 options(contrasts = rep("contr.sum", 2))

 nloglik <- function(phat,phat0) -sum(counts * log(c(t(phat)) * (1-alpha) + alpha * c(t(phat0))))

 max.iter <- control$max.iter
 objs <- numeric(max.iter)

 glm0 <- glm(counts ~  x + y, family = poisson, data = design)

 const <- sum(exp(predict(glm0) - coef(glm0)[1]))
 lambda0 <- log(n/const)
 pihat <- exp(predict(glm0) - coef(glm0)[1] + lambda0)
 phat0 <- c(t(matrix(pihat/n, byrow = TRUE, nr = K, ncol = L)))

 glm1 <- glm(counts ~  x * y, family = poisson, data = design)

 const <- sum(exp(predict(glm1) - coef(glm1)[1]))
 lambda0 <- log(n/const)
 pihat <- exp(predict(glm1) - coef(glm1)[1] + lambda0)
 phat <- c(t(matrix(pihat/n, byrow = TRUE, nr = K, ncol = L)))

 iter <- 1
 objs[iter] <- nloglik(phat, phat0)
 tol <- control$tol

 ## EM-algorithm
 while(iter < max.iter){

  # update weights
  weights_num <- phat * (1-alpha)
  weights <- weights_num / (weights_num + alpha * phat0)

  # update model parameters
  countsw <- counts * (1- weights)
  nw <- sum(countsw)
  glm0 <- glm(countsw ~  x + y, family = quasipoisson, data = design)

  const <- sum(exp(predict(glm0) - coef(glm0)[1]))
  lambda0 <- log(nw/const)
  pihat <- exp(predict(glm0) - coef(glm0)[1] + lambda0)
  phat0 <- c(t(matrix(pihat/n, byrow = TRUE, nr = K, ncol = L))) * (n/nw)
  countsw <- counts * weights
  nw <- sum(countsw)

  glm1 <- glm(countsw ~  x * y, family = quasipoisson, data = design)

  const <- sum(exp(predict(glm1) - coef(glm1)[1]))
  lambda0 <- log(nw/const)
  pihat <- exp(predict(glm1) - coef(glm1)[1] + lambda0)
  phat <- c(t(matrix(pihat/n, byrow = TRUE, nr = K, ncol = L))) * (n/nw)

  # update objective
  iter <- iter + 1
  objs[iter] <- nloglik(phat, phat0)

  if(objs[iter-1] - objs[iter] < tol)
   break

 }

 objs <- objs[1:iter]

 # calculation of covariance matrix
 ncells <- K*L
 ones <- rep(1, ncells)
 hess0 <- counts*(1-alpha)^2/((1-alpha)*phat + alpha*phat0)^2
 hess <- -(outer(ones/ncells, hess0) + outer(hess0, ones/ncells)) + mean(hess0)/ncells
 diag(hess) <- diag(hess) + hess0
 svdhess <- svd(hess)
 rootinv <- scale(svdhess$u[,1:(ncells-1)], center = FALSE, scale = sqrt(svdhess$d[1:(ncells-1)]))
 vcov_phat <- tcrossprod(rootinv)

 phat <- matrix(phat, nrow = K, ncol = L, byrow = TRUE,
                dimnames = list(unlist(attr(ctab, "row.vars")),
                                unlist(attr(ctab, "col.vars"))))
 names(dimnames(phat)) <- all.vars(formula)[2:1]
 phat0 <- matrix(phat0, nrow = K, ncol = L, byrow = TRUE,
                 dimnames = list(unlist(attr(ctab, "row.vars")),
                                 unlist(attr(ctab, "col.vars"))))
 names(dimnames(phat0)) <- all.vars(formula)[2:1]

 x <- expand.grid(names(colnames(phat)), names(rownames(phat)))
 rownames(vcov_phat) <- paste0("(", x$Var2, ", ", x$Var1, ")")
 colnames(vcov_phat) <- paste0("(", x$Var2, ", ", x$Var1, ")")

 # return output
 x <- list(objs = objs, phat = phat, phat0 = phat0, vcov_phat = vcov_phat,
           ftable = ftable(phat*n), m.rate = m.rate)
 x$call <- match.call()

 class(x) <- "ctablemixture"
 x
}
