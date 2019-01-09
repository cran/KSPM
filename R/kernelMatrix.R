#' @title Kernel matrix
#'
#' @description These functions transform a \eqn{n \times p}{n x q} matrix into a \eqn{n \times n}{n x n} kernel matrix.
#'
#'
#' @param Z a \eqn{n \times p}{n x q} matrix
#' @param whichkernel kernel function
#' @param gamma,rho,d kernel hyperparameters (see details)
#'
#' @details Given a \eqn{n \times p}{n x p} matrix, this function returns a \eqn{n \times n}{n x n} matrix where each cell represents the similarity between two samples defined by two \eqn{p-}dimensional vectors \eqn{x} and \eqn{y},
#' \itemize{
#' \item the Gaussian kernel is defined as \eqn{k(x,y) = exp\left(-\frac{\parallel x-y \parallel^2}{\rho}\right)}{k(x,y) = exp(-||x-y||^2 / rho)} where \eqn{\parallel x-y \parallel}{||x-y||} is the Euclidean distance between \eqn{x} and \eqn{y} and \eqn{\rho > 0}{rho > 0} is the bandwidth of the kernel,
#' \item the linear kernel is defined as \eqn{k(x,y) = x^Ty}{k(x,y) = t(x).y},
#' \item the polynomial kernel is defined as \eqn{k(x,y) = (\rho x^Ty + \gamma)^d}{k(x,y) = (rho.t(x).y + gamma)^d} with \eqn{\rho > 0}{rho > 0}, \eqn{d} is the polynomial order. Of note, a linear kernel is a polynomial kernel with \eqn{\rho = d = 1}{rho = d = 1} and \eqn{\gamma = 0}{gamma = 0},
#' \item the sigmoid kernel is defined as \eqn{k(x,y) = tanh(\rho x^Ty + \gamma)}{k(x,y) = tanh(rho.t(x).y + gamma)} which is similar to the sigmoid function in logistic regression,
#' \item the inverse quadratic function defined as \eqn{k(x,y) = \frac{1}{\sqrt{\parallel x-y \parallel^2 + \gamma}}}{k(x,y) = 1 / sqrt( ||x-y||^2 + gamma)} with \eqn{\gamma > 0}{gamma > 0},
#' \item the equality kernel defined as \eqn{k(x,y) = \left\lbrace \begin{array}{ll} 1 & if  x = y \\ 0 & otherwise \end{array}\right.}{k(x,y) = 1 if x = y, 0 otherwise}.
#' }
#'
#' @return A \eqn{n \times n}{n x n} matrix.
#'
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kernel.gaussian}, \link{kernel.linear}, \link{kernel.polynomial}, \link{kernel.equality}, \link{kernel.sigmoid}, \link{kernel.inverse.quadratic}.
#'
#'
#' @export



kernel.matrix <- function(Z, whichkernel, rho = NULL, gamma = NULL, d = NULL)
{

  # gaussian kernel
  if ( whichkernel == "gaussian" )
  {
    if (!is.null(rho)) {
      K <- kernel.gaussian(x = Z, rho = rho)
    } else {
      n <- dim(Z)[1]
      names <- rownames(Z)
      K <- matrix(NA, nrow = n, ncol = n)
      rownames(K) <- names
      colnames(K) <- names
    }
  }

  # linear kernel
  if ( whichkernel == "linear" )
  {
    K <- kernel.linear(x = Z)
  }

  # polynomial kernel
  if ( whichkernel == "polynomial" )
  {
    if (!is.null(rho) & !is.null(gamma)) {
      K <- kernel.polynomial(x = Z, rho = rho, gamma = gamma, d = d)
    } else {
      n <- dim(Z)[1]
      names <- rownames(Z)
      K <- matrix(NA, nrow = n, ncol = n)
      rownames(K) <- names
      colnames(K) <- names
    }
  }

  # sigmoid kernel
  if ( whichkernel == "sigmoid" )
  {
    if (!is.null(rho) & !is.null(gamma)) {
      K <- kernel.sigmoid(x = Z, rho = rho, gamma = gamma)
    } else {
      n <- dim(Z)[1]
      names <- rownames(Z)
      K <- matrix(NA, nrow = n, ncol = n)
      rownames(K) <- names
      colnames(K) <- names
    }
  }

  # inverse quadratic kernel
  if ( whichkernel == "inverse.quadratic" )
  {
    if (!is.null(rho) & !is.null(gamma)) {
      K <- kernel.inverse.quadratic(x = Z, gamma = gamma)
    } else {
      n <- dim(Z)[1]
      names <- rownames(Z)
      K <- matrix(NA, nrow = n, ncol = n)
      rownames(K) <- names
      colnames(K) <- names
    }
  }

  # equality kernel
  if ( whichkernel == "equality" )
  {
    K <- kernel.equality(x = Z)
  }

  return(K)
}
