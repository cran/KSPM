#' @title Summarizing Kernel Semi parametric Model Fits
#'
#' @title Summary method for kernel semi parametric model
#'
#' @description summary method for an object of class "kspm"
#'
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param kernel.test vector of characters indicating for which kernel a test should be performed. Default is \code{"all"}. If \code{"none"}, no test will be performed.
#' @param global.test logical, if \code{TRUE}, a global test for kernel part is computed.
#' @param ... further arguments passed to or from other methods.
#'
#' @details the description of the model, including coefficients for the linear part and if asked for, test(s) of variance components associated with kernel part.
#'
#'
#' @return Computes and returns the followimg summary statistics of the fitted kernel semi parametric model given in object
#'   \item{residuals}{residuals}
#'   \item{coefficients}{a \eqn{p \times 4}{p x 4} matrix with columns for the estimated coefficient, its standard error, t statistic and corresponding (two sided) p value for the linear part of the model.}
#'   \item{sigma}{the square root of the estimated variance of the random error \eqn{\sigma^2 = \frac{RSS}{edf} }{sigma^2 = RSS / edf} where \eqn{RSS}{RSS} is the residual sum of squares and \eqn{edf}{edf} is the effective degree of freedom.}
#'   \item{edf}{effective degrees of freedom}
#'   \item{r.squared}{\eqn{R^2}{R^2}, the fraction of variance explained by the model, \eqn{1 - \frac{\sum e_i^2}{\sum(y_i - y^{\ast})^2}}{1 - sum(e_i^2) / sum((y_i - y^star)^2)} where \eqn{y^{\ast}}{y^star} is the mean of \eqn{y_i}{y_i} if there is an intercept and zero otherwise.}
#'   \item{adj.r.squared}{the above \eqn{R^2}{R^2} statistics, adjusted, penalizing for higher \eqn{p}.}
#'   \item{score.test}{a \eqn{q \times 3}{q x 3} matrix with colums for the estimated lambda, tau and p value for the q kernels for which a test should be performed.}
#'   \item{global.p.value}{p value from the score test for the global model.}
#'   \item{sample.size}{sample size (all: global sample size, inc: complete data sample size).}
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @references
#'
#' Liu, D., Lin, X., and Ghosh, D. (2007). Semiparametric regression of multidimensional genetic pathway data: least squares kernel machines and linear mixed models. Biometrics, 63(4), 1079:1088.
#'
#' Schweiger, Regev, et al. "RL SKAT: an exact and efficient score test for heritability and set tests." Genetics (2017): genetics 300395.
#'
#' Li, Shaoyu, and Yuehua Cui. "Gene centric gene gene interaction: A model based kernel machine method." The Annals of Applied Statistics 6.3 (2012): 1134:1161.
#'
#' @examples
#' x <- 1:15
#' z1 <- runif(15, 1, 6)
#' z2 <- rnorm(15, 1, 2)
#' y <- 3*x + (z1 + z2)^2 + rnorm(15, 0, 2)
#' fit <- kspm(y, linear = ~ x, kernel = ~ Kernel(~ z1 + z2,
#' kernel.function = "polynomial", d= 2, rho = 1, gamma = 0))
#' summary(fit)
#'
#' @seealso \link{kspm} for fitting model, \link{predict.kspm} for predictions, \link{plot.kspm} for diagnostics
#'
#' @importFrom stats pt
#'
#' @rdname summary.kspm
#' @export summary.kspm
#' @export





summary.kspm <- function(object, kernel.test = "all", global.test = FALSE, ...) {

  ################################################################################
  # CALL
  ################################################################################

  call <- object$call

  ################################################################################
  # SAMPLE SIZE
  ################################################################################

  # Total sample size
  n.total <- object$n.total
  # Sample size used by the model (NA excluded)
  n <- object$n

  ################################################################################
  # LINEAR PART
  ################################################################################

  if (dim(object$X)[2] != 0) {
    # Estimates
    beta.est <- object$linear.coefficients

    # Variance-covariance matrix of estimates
    beta.vcov <- object$sigma ^ 2 * solve(t(object$X) %*% object$L %*% object$X) %*% t(object$X) %*% object$L %*% object$L %*% object$X %*% solve(t(object$X) %*% object$L %*% object$X)
    # Standard error as diagonal of the Variance-covariance matrix
    beta.se <- sqrt(diag(beta.vcov))

    # T-values
    t.value <- beta.est/beta.se

    # P-values
    p.value <- 2 * pt(abs(t.value), object$edf, lower.tail = FALSE)

    # T-Table
    coefficients <- matrix(cbind(beta.est, beta.se, t.value, p.value), nrow = length(beta.est))
    colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    rownames(coefficients) <- c("(Intercept)", colnames(object$X)[-1])
  }

  ################################################################################
  # KERNEL PART
  ################################################################################

  if (!is.null(object$K)) {
    # kernel tested should be in the model
    if (sum(!(kernel.test %in% c("all", "none", names(object$lambda)))) > 0) {
      stop("kernel.test includes a kernel not present in the model")
    }

    if (!("none" %in% kernel.test)) {
      # Test for model with one kernel
      if (length(object$lambda) == 1) {
        kernel.test <- names(object$lambda)
        tau <- (object$sigma ^ 2) / object$lambda
        res_test1 <- test.1.kernel(object)
        pTau <- res_test1$p.value
        lambda <- object$lambda
        daviesParameters <- list(Ker1 = list(q = res_test1$q, pev = res_test1$pev))
      }
      # Test for model with 2 kernels and more
      if (length(object$lambda) > 1) {
        tau <- c()
        pTau <- c()
        lambda <- c()
        daviesParameters <- list()
        # If user did not indicate which kernel should be tested, all are tested
        if ("all" %in% kernel.test) {
          kernel.test <- names(object$lambda)
        }
        for (k in 1:length(kernel.test)) {
          tau[k] <- (object$sigma ^ 2) / object$lambda[kernel.test[k]]
          res_testk <- test.k.kernel(object, kernel.test[k])
          pTau[k] <- res_testk$p.value
          lambda[k] <- object$lambda[kernel.test[k]]
          daviesParameters[[k]] <- list(q = res_testk$q, pev = res_testk$pev)

        }
        names(daviesParameters) <- kernel.test
      }
      # Table
      resKernel <- matrix(cbind(lambda, tau, pTau), nrow = length(kernel.test))
      colnames(resKernel) <- c("lambda", "tau", "p-value")
      rownames(resKernel) <- kernel.test
    } else {# If user indicated none, no test for kernel is displayed
      resKernel <- NULL
      daviesParameters <- NULL
    }

    # Global test
    if (global.test) {
      res_test_global <- test.global.kernel(object)
      pGlobal <- res_test_global$p.value
      daviesParametersGlobal <- list(q = res_test_global$q, pev = res_test_global$pev)
    } else {
      pGlobal <- NULL
      daviesParametersGlobal <- NULL
    }
  } else {
    resKernel <- NULL
    pGlobal <- NULL
    daviesParameters <- NULL
    daviesParametersGlobal <- NULL
  }

  ################################################################################
  # R2/R2adj
  ################################################################################

  # R square
  R2 <- (sum( (object$Y - mean(object$Y)) ^ 2) - sum(object$residuals ^ 2)) / sum( (object$Y - mean(object$Y)) ^ 2)
  # Adjusted R square
  R2adj <- ((sum( (object$Y - mean(object$Y)) ^ 2))/(object$n - 1) - sum(object$residuals ^ 2)/object$edf) / (sum( (object$Y - mean(object$Y)) ^ 2)/(object$n - 1))

  ################################################################################
  # OUT
  ################################################################################

  out <- list(call = call, residuals = object$residuals, coefficients = coefficients, sigma = object$sigma, edf = object$edf, R2 = R2, R2adj = R2adj, score.test = resKernel, global.p.value = pGlobal, sample.size = list(all = n.total, inc = n), daviesInfo = daviesParameters, daviesInfoGlobal = daviesParametersGlobal)
  class(out) <- c("summary.kspm")
  return(out)
}





