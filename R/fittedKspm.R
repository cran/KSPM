#' @title Extract Model Fitted values
#'
#' @description Returns fitted values for a model of class "kspm".
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param ... additional optional argument (currently unused).
#'
#'
#'
#' @return The vector of fitted values.
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @references Liu, D., Lin, X., and Ghosh, D. (2007). Semiparametric regression of multidimensional genetic pathway data: least squares kernel machines and linear mixed models. Biometrics, 63(4), 1079:1088.
#'
#' @seealso \link{kspm} for fitting model, \link{residuals.kspm}, \link{coef.kspm}, \link{nobs.kspm}.
#'
#' @examples
#' x <- 1:15
#' z <- runif(15, 1, 6)
#' y <- 3*x + z^2 + rnorm(15, 0, 2)
#' fit <- kspm(y, linear = ~ x, kernel = ~ Kernel(z,
#' kernel.function = "polynomial", d = 2, rho = 1, gamma = 0))
#' fitted(fit)
#'
#' @importFrom stats fitted
#'
#'
#' @rdname fitted.kspm
#' @export fitted.kspm
#' @export


fitted.kspm <- function(object, ...)
{
  return(object$fitted.values)
}

