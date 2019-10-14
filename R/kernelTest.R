#' @title Score Tests for kernel part in kernel semi parametric model
#'
#' @description Perform score tests for kernel part in kernel semi parametric model
#'
#' @name test.function
#' @rdname test.function
#' @aliases test.1.kernel
#' @aliases test.global.kernel
#' @aliases test.k.kernel
#'
#' @param object an object of class "kspm"
#' @param kernel.name vector of character listing names of kernels for which test should be performed
#'
#'
#' @return p values
#'
#' @references
#'
#' Schweiger, Regev, et al. "RL SKAT: an exact and efficient score test for heritability and set tests." Genetics (2017): genetics 300395.
#'
#' Li, Shaoyu, and Yuehua Cui. "Gene centric gene gene interaction: A model based kernel machine method." The Annals of Applied Statistics 6.3 (2012): 1134:1161.
#'
#' Oualkacha, Karim, et al. "Adjusted sequence kernel association test for rare variants controlling for cryptic and family relatedness." Genetic epidemiology 37.4 (2013): 366:376.
#'
#' Ge, Tian, et al. "A kernel machine method for detecting effects of interaction between multidimensional variable sets: An imaging genetics application." Neuroimage 109 (2015): 505:514.
#'
#'
#' @importFrom stats lm resid
#' @importFrom CompQuadForm davies
#' @importFrom expm sqrtm
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @export



################################################################################
# TEST OF KERNEL IN SINGLE KSPM
################################################################################

test.1.kernel <- function(object)
{

  #-------------------------------------------------------------------------------
  # Information from the model
  #-------------------------------------------------------------------------------

  # sample size
  n <- object$n
  # outcome
  Y <- object$Y
  # linear part
  X <- object$X
  # tested kernel
  K <- object$K[[1]]

  # if the model under H0 is null without intercept, no test
  if (dim(X)[2] == 0) {
    p.value <- NaN
    warning("test cannot be performed: no variable, no intercept in the null model")
  } else {

    #-------------------------------------------------------------------------------
    # Model under H0
    #-------------------------------------------------------------------------------

    # model
    mod <- lm(Y ~ -1 + X)

    #-------------------------------------------------------------------------------
    # Test
    #-------------------------------------------------------------------------------

    # residuals
    res <- resid(mod)
    # variance of residuals
    s2 <- sum(res ^ 2) # OU s2 <- mod$sigma^2
    # test (exact method)
    I <- diag(n)
    S <-  (I - X %*% solve(t(X) %*% X) %*% t(X))
    q <- as.numeric(res %*% K %*% res / s2)
    SKS <- S %*% K %*% S - q * S
    ee <- eigen(SKS, symmetric = T, only.values = TRUE)
    pev <- ee$values[abs(ee$values) >= 1e-10]
    davies_res <- davies(q, lambda = pev)
    k_lim <- 1
    while (davies_res$ifault == 1 & k_lim <= 5) {
      k_lim <- k_lim + 1
      davies_res <- davies(q, lambda = pev, lim = 10000*k_lim)
    }
    p.value <- davies_res$Qq
  }

  return(list(q = q, pev = pev, p.value = p.value))
}

################################################################################
# GLOBAL TEST OF KERNEL
################################################################################

#' @rdname test.function
#' @export
test.global.kernel <- function(object)
{

  #-------------------------------------------------------------------------------
  # Information from the model
  #-------------------------------------------------------------------------------

  # sample size
  n <- object$n
  # outcome
  Y <- object$Y
  # linear part
  X <- object$X
  # tested kernel
  K <- Reduce("+", object$K)

  # if the model under H0 is null without intercept, no test
  if (dim(X)[2] == 0) {
    p.value <- NaN
    warning("global test cannot be performed: no variable, no intercept in the null model")
  } else {

    #-------------------------------------------------------------------------------
    # Model under H0
    #-------------------------------------------------------------------------------

    # model
    mod <- lm(Y ~ -1 + X)

    #-------------------------------------------------------------------------------
    # Test
    #-------------------------------------------------------------------------------

    # residuals
    res <- resid(mod)
    # variance of residuals
    s2 <- sum(res ^ 2) # OU s2 <- mod$sigma^2
    # test (exact method)
    I <- diag(n)
    S <-  (I - X %*% solve(t(X) %*% X) %*% t(X))
    q <- as.numeric(res %*% K %*% res / s2)
    SKS <- S %*% K %*% S - q * S
    ee <- eigen(SKS, symmetric = T, only.values = TRUE)
    pev <- ee$values[abs(ee$values) >= 1e-10] # positive eigen values
    davies_res <- davies(q, lambda = pev)
    k_lim <- 1
    while (davies_res$ifault == 1 & k_lim <= 5) {
      k_lim <- k_lim + 1
      davies_res <- davies(q, lambda = pev, lim = 10000*k_lim)
    }
    p.value <- davies_res$Qq
  }

  return(list(q = q, pev = pev, p.value = p.value))
}



################################################################################
# TEST OF KERNEL IN MULTIPLE KSPM
################################################################################

#' @rdname test.function
#' @export
test.k.kernel <- function(object, kernel.name)
{

  #-------------------------------------------------------------------------------
  # Information from the model
  #-------------------------------------------------------------------------------

  # sample size
  n <- object$n
  # outcome
  Y <- object$Y
  # linear part
  X <- object$X
  # tested kernel
  K.H1 <- object$K[kernel.name][[1]]

  #-------------------------------------------------------------------------------
  # Model under H0
  #-------------------------------------------------------------------------------

  # kernel part under H0
  K.H0 <- object$K[!(names(object$K) %in% kernel.name)]
  # model under H0
  param <- search.parameters(Y = Y, X = X, kernelList = K.H0,  n = n, not.missing = NULL, compute.kernel = FALSE, controlKspm = kspmControl())
  # fitted values
  fitted.Kalpha <- list()
  for (i in 1:length(param$K)) {
    fitted.Kalpha[[i]] <- param$K[[i]] %*% param$alpha[[i]]
  }
  if (dim(X)[2] == 0) {
    fitted.values <- Reduce("+", fitted.Kalpha)
  } else{
    fitted.values <- X %*% param$beta + Reduce("+", fitted.Kalpha)
  }

  # residuals
  res <- Y - fitted.values

  # effective degree of freedom
  edf <- sum(diag(param$Hat))
  sigma <- sqrt((1 / (n - edf)) * sum(res ^ 2)) # DEMANDER AVIS A CELIA ET AURELIE POUR CETTE FORMULE

  #-------------------------------------------------------------------------------
  # Test
  #-------------------------------------------------------------------------------

  # sigma matrix under H0
  I <- diag(n)
  lambda_invK = list()
  for (i in 1:length(K.H0)) {
    lambda_invK[[i]] <- 1/param$lambda[i] * param$K[[i]]
  }
  sum_lambda_invK <- Reduce("+", lambda_invK)
  sigma0 <- sigma^2 * (sum_lambda_invK + I)
  sigma0_inv <- solve(sigma0)

  # if no X
  if (dim(X)[2] == 0) {
    P0 <- sigma0_inv
  } else {
    # if X
    P0 <- sigma0_inv %*% (I - X %*% solve(t(X) %*% sigma0_inv %*% X) %*% t(X) %*% sigma0_inv)
  }
  # test statistic
  q <- (1 / 2) * as.numeric(t(Y) %*% P0  %*% K.H1 %*% P0 %*% Y)
  # distribution under H0
  P0.sqrt <- sqrtm(P0) # root square of the matrix
  ee <- eigen((1 / 2) * t(P0.sqrt) %*% K.H1 %*% P0.sqrt,  only.values = TRUE)
  pev <- ee$values[abs(ee$values) >= 1e-10]
  davies_res <- davies(q, lambda = pev)
  k_lim <- 1
  while (davies_res$ifault == 1 & k_lim <= 5) {
    k_lim <- k_lim + 1
    davies_res <- davies(q, lambda = pev, lim = 10000*k_lim)
  }
  p.value <- davies_res$Qq
  return(list(q = q, pev = pev, p.value = p.value))
}
