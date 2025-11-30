local_edition(3)

test_that("glm_ele example works", {

  ########################################################################
  ## 1. Generate the Test Data
  ########################################################################

  set.seed(4)

  Afun <- function(N, lambda) {
    A <- matrix(0, nrow = N, ncol = N)
    set <- NULL

    for(i in 1:N) {
      r <- runif(1)
      if(r <= lambda) {
        set <- c(set, i)
        A[i, i] <- 1
      }
    }

    wrong <- (1:N)[-set]
    for(i in wrong) {
      valid <- unique((1:N)[-c(i, set)])
      if(length(valid) > 1) j <- sample(valid, size = 1)
      else j <- valid
      A[i, j] <- 1
      set <- c(set, j)
    }

    A
  }

  Q <- 3
  M <- c(150, 350, 500)
  N <- sum(M)
  Truebeta <- c(1, -1, 0.5, 1.1, 1.5)
  p <- length(Truebeta)
  Truelambda <- c(0.99, 0.8, 0.7)
  Estlambda <- Truelambda

  Block <- c(rep(1, M[1]), rep(2, M[2]), rep(3, M[3]))
  block <- as.factor(Block)

  contrasts(block) <- contr.sum(length(levels(block)))

  tdm <- t(cbind(1, contrasts(block)))
  mr_calc <- function(x) { c(plogis(x %*% tdm)) }

  Truegamma <- c(solve(tdm %*% t(tdm)) %*% tdm %*% qlogis(Truelambda))
  Truesigma <- 0.25

  # generate X and Y
  data <- data.frame(
    X0 = rep(1, N),
    X1 = rnorm(N, 0, 1),
    X2 = rnorm(N, 1, 1),
    block = block
  )

  X <- model.matrix(~ X1 + X2 + block, data = data)
  L <- c(X %*% matrix(Truebeta, ncol = 1))
  Y <- rnorm(N, L, Truesigma^2)

  # permute Y with linkage errors
  Ystar <- c()
  for(q in 1:Q) {
    Aq <- Afun(N = M[q], lambda = Truelambda[q])
    Ystar <- c(Ystar, c(Aq %*% Y[Block == q]))
  }

  ########################################################################
  ## 2. Fit model using glm_ele
  ########################################################################

  expect_snapshot({
    fit <- glm_ele(Ystar ~ X - 1, m.rate = 1 - Estlambda, blocks = Block)
    fit
  })

  expect_snapshot(fit$coefficients)
  expect_snapshot(fit$standard.errors)
  expect_snapshot(fit$dispersion)
  expect_snapshot(summary(fit))
  expect_snapshot(vcov(fit))
  expect_snapshot(confint(fit))
  expect_snapshot(predict(fit))

})
