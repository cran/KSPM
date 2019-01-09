#' @title Extract the number of observations from a Kernel Semi parametric Model Fit
#'
#' @description Extract the number of observations use to estimate the model coefficients. This is principally intented to be used in computing BIC (see \link{extractAIC.kspm}).
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param ... additional optional argument (currently unused).
#'
#'
#'
#' @return A single number (integer).
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kspm} for fitting model, \link{extractAIC.kspm}.
#'
#' @examples
#' x <- 1:15
#' y <- 3*x + rnorm(15, 0, 2)
#' fit <- kspm(y, kernel = ~ Kernel(x, kernel.function = "linear"))
#' nobs(fit)
#'
#' @importFrom stats nobs
#'
#' @rdname nobs.kspm
#' @export nobs.kspm
#' @export


nobs.kspm <- function(object, ...)
{
  return(object$n)
}


