#' @title Extract Model Hyper-parameter
#'
#' @description Returns hyper-parameters for a model of class "kspm".
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param ... additional optional argument (currently unused).
#'
#'
#'
#' @return A list of parameter.
#' \item{lambda}{A vector of penalisation arameters.}
#' \item{kernel}{A vector of tunning parameters.}
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
#' hypercoef(fit)
#'
#'
#' @rdname hypercoef
#' @export hypercoef
#' @export


hypercoef <- function(object, ...)
{
  lambda <- object$lambda
  tuning <- c(ifelse(is.null(object$kernel.info[[1]]$rho), NA, object$kernel.info[[1]]$rho), ifelse(is.null(object$kernel.info[[1]]$d), NA, object$kernel.info[[1]]$d), ifelse(is.null(object$kernel.info[[1]]$gamma), NA, object$kernel.info[[1]]$gamma))
  names(tuning) <- c("rho", "d", "gamma")
  if( length(object$kernel.info) > 1) {
    for (i in 2:length(object$kernel.info)) {
      tuning <- rbind(tuning, c(ifelse(is.null(object$kernel.info[[i]]$rho), NA, object$kernel.info[[i]]$rho), ifelse(is.null(object$kernel.info[[i]]$d), NA, object$kernel.info[[i]]$d), ifelse(is.null(object$kernel.info[[i]]$gamma), NA, object$kernel.info[[i]]$gamma)))
    }
    row.names(tuning) <- names(object$kernel.info)
  } else {
    tuning <- tuning[!is.na(tuning)]
  }
  out <- list(lambda = lambda, tuning = tuning)
  return(out)
}



