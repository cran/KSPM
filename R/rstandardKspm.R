#' @title Standardized residuals for Kernel Semi parametric Model Fits
#'
#' @description computes standardized residuals for an object of class "kspm".
#'
#'
#'
#' @param model an model of class "kspm", usually, a result of a call to \code{kspm}.
#' @param ... furter arguments passed to or from other methods (currently unused).
#'
#' @details Standardized residuals \eqn{t_i}{t_i} are obtained by \eqn{t_i = \frac{e_i}{\hat{\sigma} \sqrt{1 - h_{ii}}}}{t_i = e_i / hat(sigma) sqrt(1 - h_ii)} where \eqn{e_i}{e_i} is the residual, \eqn{\hat{\sigma}}{hat(sigma)} is the estimated standard deviation of the errors and \eqn{h_{ii}}{h_ii} is the leverage of subject i, i.e. the i th diagonal element of the Hat matrix.
#'
#'
#' @return a vector containing the standardized residuals.
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kspm} for fitting model, \link{residuals.kspm}, \link{cooks.distance.kspm}, \link{plot.kspm}.
#'
#' @importFrom stats rstandard
#'
#'
#' @rdname rstandard.kspm
#' @export rstandard.kspm
#' @export


rstandard.kspm <- function(model, ...)
{
  # Leverage is digonal of Hat matrix
  leverage <- diag(model$Hat)
  # Compute standardized residuals
  standardized.residuals <- model$residuals / (model$sigma * sqrt(1 - leverage))
  return(standardized.residuals)
}

