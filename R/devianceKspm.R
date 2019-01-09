#' @title Model deviance
#'
#' @description Returns the deviance of a fitted model object of class "kspm".
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}, for which the deviance is desired.
#' @param ... additional optional argument (currently unused).
#'
#' @details This function extracts deviance of a model fitted using \code{kspm} function. The returned deviance is the residual sum of square (RSS).
#'
#'
#' @return The value of the deviance extracted from the object \code{object}.
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kspm}, \link{extractAIC.kspm}
#'
#' @examples
#' x <- 1:15
#' y <- 3*x + rnorm(15, 0, 2)
#' fit <- kspm(y, kernel = ~ Kernel(x, kernel.function = "linear"))
#' deviance(fit)
#'
#' @importFrom stats deviance
#'
#' @rdname deviance.kspm
#' @export deviance.kspm
#' @export


deviance.kspm <- function(object, ...)
{
  # Compute the Residual Sum of Squares
  RSS <- (t(object$residuals) %*% (object$residuals))[1, 1]
  # Deviance is the RSS in Gaussian family models
  cat("RSS")
  return(RSS)
}



