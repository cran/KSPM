#' @title Extract Model Coefficients
#'
#' @description Returns linear and kernel coefficients for a model of class "kspm".
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param ... additional optional argument (currently unused).
#'
#'
#'
#' @return Two matrices of coefficients.
#' \item{linear}{A vector of coefficients for linear part. One row is one variable.}
#' \item{kernel}{A matrix of coefficients for linear part. One row is one subject, one column is one kernel part.}
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @references Liu, D., Lin, X., and Ghosh, D. (2007). Semiparametric regression of multidimensional genetic pathway data: least squares kernel machines and linear mixed models. Biometrics, 63(4), 1079:1088.
#'
#' @seealso \link{kspm} for fitting model.
#'
#' @examples
#' x <- 1:15
#' z1 <- runif(15, 1, 6)
#' z2 <- rnorm(15, 1, 2)
#' y <- 3*x + (z1 + z2)^2 + rnorm(15, 0, 2)
#' fit <- kspm(y, linear = ~ x, kernel = ~ Kernel(~ z1 + z2,
#' kernel.function = "polynomial", d= 2, rho = 1, gamma = 0))
#' coef(fit)
#'
#' @importFrom stats coef
#'
#' @rdname coef.kspm
#' @export coef.kspm
#' @export


coef.kspm <- function(object, ...)
{
  linear <- object$linear.coefficients
  colnames(linear) <- "estimates"
  out <- list(linear = linear, kernel = object$kernel.coefficients)
  return(out)
}



