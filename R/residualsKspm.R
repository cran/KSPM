#' @title Extract residuals from a Kernel Semi Parametric Model
#'
#' @description Returns the vector of residuals for a model fit of class "kspm".
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param ... additional optional argument (currently unused).
#'
#'
#'
#' @return A vector of residuals. The vector length is the number of observations used in model coefficients estimation (see \link{nobs.kspm}).
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kspm} for fitting model, \link{nobs.kspm}, \link{rstandard.kspm}.
#'
#' @examples
#' x <- 1:15
#' y <- 3*x + rnorm(15, 0, 2)
#' fit <- kspm(y, kernel = ~ Kernel(x, kernel.function = "linear"))
#' residuals(fit)
#'
#' @importFrom stats residuals
#'
#' @rdname residuals.kspm
#' @export residuals.kspm
#' @export




residuals.kspm <- function(object, ...)
{
  return(object$residuals)
}


