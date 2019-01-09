#' @title Confidence interavls for linear part of model parameters
#'
#' @description Computes confidence intervals for one or more parameters in the linear part of a fitted model of class "kspm".
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param parm a vector of names specifying which parameters are to be given confidence intervals. If missing, all parameters are considered.
#' @param level the confidence level required. By default 0.95.
#' @param ... additional optional argument (currently unused).
#'
#' @details For objects of class "kspm", the confidence interval is based on student distribution and effective degree of freedom of the model.
#'
#'
#' @return A matrix with column giving lower and upper confidence limits for each parameter. These are labelled as \eqn{\frac{1-level}{2}}{(1-level) / 2} and \eqn{1 - \frac{1-level}{2}}{1-(1-level)/2} in percentage.
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kspm} for fitting model, \link{summary.kspm}.
#'
#' @examples
#' x <- 1:15
#' z1 <- runif(15, 1, 6)
#' z2 <- rnorm(15, 1, 2)
#' y <- 3*x + (z1 + z2)^2 + rnorm(15, 0, 2)
#' fit <- kspm(y, linear = ~ x, kernel = ~ Kernel(~ z1 + z2,
#' kernel.function = "polynomial", d= 2, rho = 1, gamma = 0))
#' confint(fit)
#'
#' @importFrom stats confint
#'
#'
#' @rdname confint.kspm
#' @export confint.kspm
#' @export


confint.kspm <- function(object, parm = NULL, level = 0.95, ...)
{
  if (is.null(parm)) {
    parm <- rownames(object$linear.coefficients)
  } else {
    if (sum(parm %in% rownames(object$linear.coefficients)) > 0) {
      stop(paste(parm[!(parm %in% rownames(object$linear.coefficients))][1], "is not a variable of linear part of the model"))
    }
  }
  # Estimates
  estimates <- as.matrix(object$linear.coefficients[parm, ])
  # Variance-covariance matrix of estimates
  beta.vcov <- object$sigma ^ 2 * solve(t(object$X) %*% object$L %*% object$X) %*% t(object$X) %*% object$L %*% object$L %*% object$X %*% solve(t(object$X) %*% object$L %*% object$X)
  # Standard error as diagonal of the Variance-covariance matrix
  estimates.SE <- cbind(estimates, as.matrix(sqrt(diag(beta.vcov))))
  # Colnames
  colnames(estimates.SE) <- c("Estimates", "SE")
  # Compute quantile from Student distribution at level "level" and "edf" effective degree of freedom
  St.quantile <- qt(p = (1-level)/2, df = object$edf)
  # Compute confidence interval of estimates
  confint.lower <- as.matrix(estimates.SE[, "Estimates"] + St.quantile * estimates.SE[, "SE"])
  confint.upper <- as.matrix(estimates.SE[, "Estimates"] - St.quantile * estimates.SE[, "SE"])
  # Return matrix of confidence interval
  out <- cbind(confint.lower, confint.upper)
  colnames(out) <- c(paste(100*(1-level)/2, "%"), paste(100 - 100*(1-level)/2, "%"))
  return(out)
}


