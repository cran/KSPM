#' @title Summarizing Kernel Semi parametric Model Fits with flexible parameters for Davies' approximation method
#'
#' @title Summary of kernel semi parametric model after adjusting parameters for Davies' approximation method
#'
#' @description for flexibility in summary method for an object of class "summary.kspm"
#'
#'
#'
#'
#' @param object an object of class "summary.kspm", usually, a result of a call to \code{summary.kspm}.
#' @param method method to approximate the chi square distribution in p-value computation, default is 'davies', another possibility is 'imhof'.
#' @param acc,lim see davies and imhof functions in CompQuadForm package.
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
#' summary.fit <- summary(fit)
#' flexible.summary(summary.fit, acc = 0.000001, lim = 1000)
#'
#' @seealso \link{kspm} for fitting model, \link{predict.kspm} for predictions, \link{plot.kspm} for diagnostics
#'
#' @importFrom CompQuadForm davies imhof
#'
#' @rdname flexible.summary
#' @export flexible.summary
#' @export





flexible.summary <- function(object, method = "davies", acc = 0.000001, lim = 10000) {

  resKernel <- object$score.test
  pGlobal <- object$pGlobal

  if (is.null(resKernel) & is.null(pGlobal)) {
    warning("No kernel part tested to apply flexible.summary")
  } else {
    if (!is.null(resKernel)) {
      for (i in row.names(resKernel)) {
        if (method == "davies") {
          davies_res <- davies(object$daviesInfo[[i]]$q, lambda = object$daviesInfo[[i]]$pev, acc = acc, lim = lim)
          resKernel[i, "p-value"] <- davies_res$Qq
        } else if (method == "imhof") {
          imhof_res <- imhof(object$daviesInfo[[i]]$q, lambda = object$daviesInfo[[i]]$pev, epsabs = acc, limit = lim)
          resKernel[i, "p-value"] <- imhof_res$Qq
        }

      }
    }
    if (!is.null(pGlobal)) {
      if (method == "davies") {
        davies_res_global <- davies(object$daviesInfoGlobal$q, lambda = object$daviesInfoGlobal$pev, acc = acc, lim = lim)
        pGlobal <- davies_res_global$Qq
      } else if (method == "imhof") {
        imhof_res <- imhof(object$daviesInfoGlobal[[i]]$q, lambda = object$daviesInfoGlobal[[i]]$pev, epsabs = acc, limit = lim)
        pGlobal <- imhof_res$Qq
      }

    }
  }




  ################################################################################
  # OUT
  ################################################################################

  out <- list(call = object$call, residuals = object$residuals, coefficients = object$coefficients, sigma = object$sigma, edf = object$edf, R2 = object$R2, R2adj = object$R2adj, score.test = resKernel, global.p.value = pGlobal, sample.size = object$sample.size, daviesInfo = object$daviesParameters, daviesInfoGlobal = object$daviesParametersGlobal)
  class(out) <- c("summary.kspm")
  return(out)
}





