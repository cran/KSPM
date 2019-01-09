#' @title compute Kernel Semi Parametric model parameters
#'
#' @description internal function to compute model parameters
#'
#' @param X X matrix
#' @param Y response matrix
#' @param kernelList list of kernels
#' @param free.parameters free parameters
#' @param n number of samples
#' @param not.missing number of non missing samples
#' @param compute.kernel boolean indicating if kernel should be computed
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @importFrom stats sd
#'
#' @export





get.parameters <- function(X = NULL, Y = NULL, kernelList = NULL, free.parameters = NULL, n = NULL, not.missing = NULL, compute.kernel = NULL)
{


  ################################################################################
  # MODIFY KERNEL MATRIX IF NEEDED
  ################################################################################

  # total number of tuning parameters
  free.parameters.count <- 1

  # if compute.kernel = TRUE
  if (compute.kernel) {
    # replace rho and gamma by new values
    for (i in 1:length(kernelList$listKernel)) {
      if ("rho" %in% kernelList$listKernel[[i]]$free.parameters) {
        kernelList$listKernel[[i]]$rho <- free.parameters[free.parameters.count]
        free.parameters.count <- free.parameters.count + 1
      }
      if ("gamma" %in% kernelList$listKernel[[i]]$free.parameters) {
        kernelList$listKernel[[i]]$gamma <- free.parameters[free.parameters.count]
        free.parameters.count <- free.parameters.count + 1
      }
    }
    # compute new kernel list
    kernelList <- computes.KernelALL(kernelList, not.missing = not.missing)

    # keep K
    K.all <- list()
    for (i in 1:length(kernelList$listKernel)) {
      K.all[[i]] <- kernelList$listKernel[[i]]$K
    }
    K <- K.all[kernelList$inclKernel]
  } else {# if compute.kernel = FALSE
    # K.all <- list()
    # for (i in 1:length(kernelList$listKernel)) {
    #   K.all[[i]] <- kernelList$listKernel[[i]]$K
    # }
    # K <- K.all[kernelList$inclKernel]
    K <- kernelList

  }
  # replace lambda by new values
  lambda <- free.parameters[free.parameters.count:length(free.parameters)]
  # Identity matrix
  I <- diag(1, ncol = n, nrow = n)


  ################################################################################
  # COMPUTE RECURSIVELY G AND M MATRIX FOR EACH KERNEL
  ################################################################################

  # final G matrix
  G <- list()
  # final M matrix
  M <- list()

  # final K %*% G-1 %*% M matrix
  KGinvM <- list()
  # final G-1 %*% M matrix
  GinvM <- list()

  # Compute final G and M for each kernel
  if (length(K) > 1) {

    # T1<-Sys.time()
    for (i in 1:length(K)) {
      # T1<-Sys.time()
      K_ord <- K[-i]
      lambda_ord <- lambda[-i]

      # G matrix (iterations for kernel i)
      g <- list()

      # M MATRIX (iterations for kernel i)
      m <- list()

      # Initialization
      g[[1]] <- I
      m[[1]] <- I

      # Compute G and H at each (intermediate) iteration for kernel i
      for (j in 1:length(K_ord)) {
        g[[j + 1]] <- lambda_ord[j] * I + m[[j]] %*% K_ord[[j]]
        g_j_inv <- solve(g[[j + 1]])
        m[[j + 1]] <- (I - m[[j]] %*% K_ord[[j]] %*% g_j_inv) %*% m[[j]]
      }

      # Compute final G for kernel i
      G[[i]] <- lambda[i] * I + m[[length(K)]] %*% K[[i]]

      # Compute final M for kernel i
      g_i_inv <- solve(G[[i]])
      M[[i]] <- (I - m[[length(K)]] %*% K[[i]] %*% g_i_inv) %*% m[[length(K)]]

      # Compute final K %*% G-1 %*% H for kernel i
      KGinvM[[i]] <- K[[i]] %*% g_i_inv %*% m[[length(K)]]

      # Compute final G-1 %*% H for kernel i
      GinvM[[i]] <- g_i_inv %*% m[[length(K)]]
    }

  } else {

    GinvM[[1]] <- solve(lambda * I + K[[1]])
    KGinvM[[1]] <- K[[1]] %*% GinvM[[1]]
  }


  ################################################################################
  # COMPUTE LINEAR COEFFICIENTS (BETA)
  ################################################################################

  L <- I - Reduce("+", KGinvM)
  if (dim(X)[2] == 0) {
    XLX_inv <- NULL
    beta <- NULL
  } else {
    XLX_inv <- solve(t(X) %*% L %*% X)
    beta <- XLX_inv %*% t(X) %*%  L %*% Y
  }


  ################################################################################
  # COMPUTE KERNEL COEFFICIENTS (ALPHA)
  ################################################################################

  alpha <- list()
  for (i in 1:length(K)) {
    if (dim(X)[2] == 0) {
      alpha[[i]] <- GinvM[[i]] %*% Y
    } else {
      alpha[[i]] <- GinvM[[i]] %*% (Y - X %*% beta)
    }
  }


  ################################################################################
  # COMPUTE HAT MATRIX AND LOOE
  ################################################################################

  Kalpha <- list()
  for (i in 1:length(K)) {
    Kalpha[[i]] <- K[[i]] %*% alpha[[i]]
  }

  # residuals
  if (dim(X)[2] == 0) {
    res <- Y - Reduce("+", Kalpha)
  } else {
    res <- Y - X %*% beta - Reduce("+", Kalpha)
  }

  # hat matrix
  sum_KGinvM <- Reduce("+", KGinvM)
  if (dim(X)[2] == 0) {
    Hat <- sum_KGinvM
  } else {
    Hat <- sum_KGinvM %*% (I - X %*% XLX_inv %*% t(X) %*%  L) + X %*% XLX_inv %*% t(X) %*%  L
  }
  # looe
  looe <- (res / (1 - diag(Hat)))^2
  looe.mean <- mean(looe)
  looe.sem <- sd(looe) / sqrt(n)


  # return parameters
  return(list(looe.mean = looe.mean, looe.sem = looe.sem, beta = beta, alpha = alpha, K = K, Hat = Hat, L = L, XLX_inv = XLX_inv, GinvM = GinvM,lambda = lambda, kernelList = kernelList))
}



