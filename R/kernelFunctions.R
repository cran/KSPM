#' @title Kernel Functions
#'
#' @description These functions transform a \eqn{n \times p}{n x p} matrix into a \eqn{n \times n}{n x n} kernel matrix.
#'
#' @name kernel.function
#' @rdname kernel.function
#' @aliases kernel.gaussian
#' @aliases kernel.linear
#' @aliases kernel.polynomial
#' @aliases kernel.equality
#' @aliases kernel.sigmoid
#' @aliases kernel.inverse.quadratic
#'
#' @param x a \eqn{n \times p}{n x p} matrix
#' @param gamma,rho,d kernel hyperparameters (see details)
#'
#' @details Given two \eqn{p-}dimensional vectors \eqn{x} and \eqn{y},
#' \itemize{
#' \item the Gaussian kernel is defined as \eqn{k(x,y) = exp\left(-\frac{\parallel x-y \parallel^2}{\rho}\right)}{k(x,y) = exp(-||x-y||^2 / rho)} where \eqn{\parallel x-y \parallel}{||x-y||} is the Euclidean distance between \eqn{x} and \eqn{y} and \eqn{\rho > 0}{rho > 0} is the bandwidth of the kernel,
#' \item the linear kernel is defined as \eqn{k(x,y) = x^Ty}{k(x,y) = t(x).y},
#' \item the polynomial kernel is defined as \eqn{k(x,y) = (\rho x^Ty + \gamma)^d}{k(x,y) = (rho.t(x).y + gamma)^d} with \eqn{\rho > 0}{rho > 0}, \eqn{d} is the polynomial order. Of note, a linear kernel is a polynomial kernel with \eqn{\rho = d = 1}{rho = d = 1} and \eqn{\gamma = 0}{gamma = 0},
#' \item the sigmoid kernel is defined as \eqn{k(x,y) = tanh(\rho x^Ty + \gamma)}{k(x,y) = tanh(rho.t(x).y + gamma)} which is similar to the sigmoid function in logistic regression,
#' \item the inverse quadratic function defined as \eqn{k(x,y) = \frac{1}{\sqrt{\parallel x-y \parallel^2 + \gamma}}}{k(x,y) = 1 / sqrt( ||x-y||^2 + gamma)} with \eqn{\gamma > 0}{gamma > 0},
#' \item the equality kernel defined as \eqn{k(x,y) = \left\lbrace \begin{array}{ll} 1 & if  x = y \\ 0 & otherwise \end{array}\right.}{k(x,y) = 1 if x = y, 0 otherwise}.
#' }
#' Of note, Gaussian, inverse quadratic and equality kernels are measures of similarity resulting to a matrix containing 1 along the diagonal.
#'
#' @return A \eqn{n \times n}{n x n} matrix.
#'
#' @importFrom stats na.omit dist
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @references Liu, D., Lin, X., and Ghosh, D. (2007). Semiparametric regression of multidimensional genetic pathway data: least squares kernel machines and linear mixed models. Biometrics, 63(4), 1079:1088.
#'
#' @export


################################################################################
# GAUSSIAN KERNEL
################################################################################

kernel.gaussian <- function(x, rho = ncol(x))
{
  missing <- attr(na.omit(x), "na.action")
  out <- exp(-1*as.matrix(dist(x) ^ 2)/rho)
  out[, missing] <- NA
  out[missing, ] <- NA
  return(out)
}


################################################################################
# LINEAR KERNEL
################################################################################

#' @rdname kernel.function
#' @export
kernel.linear <- function(x)
{
  return(x %*% t(x))
}

################################################################################
# POLYNOMIAL KERNEL
################################################################################

#' @rdname kernel.function
#' @export
kernel.polynomial <- function(x, rho = 1, gamma = 0, d = 1)
{
  return((rho * x %*% t(x) + gamma) ^ d)
}

################################################################################
# SIGMOID KERNEL
################################################################################

#' @rdname kernel.function
#' @export
kernel.sigmoid <- function(x, rho = 1, gamma = 1)
{
  return(tanh(rho * x %*% t(x) + gamma))
}

################################################################################
# INVERSE QUADRATIC KERNEL
################################################################################

#' @rdname kernel.function
#' @export
kernel.inverse.quadratic <- function(x, gamma = 1)
{
  missing <- attr(na.omit(x), "na.action")
  out <- 1 / sqrt(as.matrix(dist(x) ^ 2) + gamma)
  out[, missing] <- NA
  out[missing, ] <- NA
  return(out)
}

################################################################################
# EQUALITY KERNEL
################################################################################

#' @rdname kernel.function
#' @export
kernel.equality <- function(x)
{
  missing <- attr(na.omit(x), "na.action")
  out <- ifelse(as.matrix(dist(x) ^ 2) == 0, 1, 0)
  out[, missing] <- NA
  out[missing, ] <- NA
  return(out)
}

