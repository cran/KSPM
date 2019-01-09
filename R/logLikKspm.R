#' @title Log Likelihood of a kspm Object
#'
#' @description Returns the Log Likelihood value of the kernel semi parametric model represented by \code{obect} evaluated at the estimated coefficients.
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \link{kspm}.
#' @param ... additional optional argument (currently unused).
#'
#' @details The function returns the Log Likelihood computed as follow:  \eqn{logLik = -\frac{1}{2} RSS}{-(1/2)*RSS} where \eqn{RSS}{RSS} is the residual sum of squares.
#'
#'
#' @return logLik of kspm fit
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @references Liu, D., Lin, X., and Ghosh, D. (2007). Semiparametric regression of multidimensional genetic pathway data: least squares kernel machines and linear mixed models. Biometrics, 63(4), 1079:1088.
#'
#' @seealso \link{kspm}, \link{extractAIC.kspm}, \link{deviance.kspm}
#'
#' @examples
#' x <- 1:15
#' y <- 3*x + rnorm(15, 0, 2)
#' fit <- kspm(y, kernel = ~ Kernel(x, kernel.function = "linear"))
#' logLik(fit)
#'
#' @importFrom stats logLik
#'
#' @rdname logLik.kspm
#' @export logLik.kspm
#' @export


logLik.kspm <- function(object, ...)
{
  # Compute the Residual Sum of Squares
  RSS <- (t(object$residuals) %*% (object$residuals))[1, 1]
  # log Lik is -(1/2)*RSS in Gaussian family models
  out <- -(1/2)*RSS
  # some attributes of the logLik class
  attr(out, "nobs") <- sum(object$n)
  attr(out, "df") <- NA # because it is edf and not df
  # class
  class(out) <- "logLik"
  return(out)
}


