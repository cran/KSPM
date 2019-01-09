#' @title Extract residuals standard deviation
#'
#' @description Returns the residuals standard deviation (sigma) for object of class "kspm".
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param ... additional optional argument (currently unused).
#'
#' @details The value returned by the method is \eqn{\sqrt{\frac{RSS}{edf}}}{sqrt(RSS / edf)} where \eqn{RSS}{RSS} is the residual sum of squares and \eqn{edf}{edf} is the effective degree of freedom.
#'
#'
#' @return typically a number, the estimated standard deviation of the errors ("residual standard deviation")
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kspm} for fitting model, \link{summary.kspm}, \link{residuals.kspm}, \link{nobs.kspm}, \link{deviance.kspm}.
#'
#'
#' @importFrom stats sigma
#'
#' @rdname sigma.kspm
#' @export sigma.kspm
#' @export


sigma.kspm <- function(object, ...)
{
  return(object$sigma)
}

