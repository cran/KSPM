#' @title Cook's distance for a Kernel Semi Parametric Model Fit
#'
#' @description Computes the Cook's distance method for an object of class "kspm".
#'
#'
#'
#' @param model an model of class "kspm", usually, a result of a call to \code{kspm}.
#' @param ... furter arguments passed to or from other methods (currently unused).
#'
#' @details Cook's distance values (\eqn{C_i}{C_i}) are computed as follows: \eqn{C_i = \frac{e_i^2 h_{ii}}{\hat{\sigma}^2 tr(H) (1-h_{ii})^2}}{C_i = (e_i^2  h_ii) / (hat(sigma^2)  tr(H)  (1-h_ii)^2)} where e_i is the residual of subject i, h_ii is the i th diagonal element of Hat matrix H corresponding to the leverage associated with subject i and tr(H) is the trace of the Hat matrix H.
#'
#'
#' @return A vector containing Cook's distance values.
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kspm} for fitting model, \link{residuals.kspm}, \link{rstandard.kspm}, \link{plot.kspm}.
#'
#' @importFrom stats cooks.distance
#'
#' @rdname cooks.distance.kspm
#' @export cooks.distance.kspm
#' @export


cooks.distance.kspm <- function(model, ...)
{
  # Leverage is digonal of Hat matrix
  leverage <- diag(model$Hat)
  # Compute Cook's distance
  cooks.distance.values <- model$residuals ^ 2 * leverage / (model$sigma ^ 2 * sum(leverage) * (1 - leverage) ^ 2)
  return(cooks.distance.values)
}


