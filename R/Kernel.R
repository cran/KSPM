#' @title Create a Kernel Object
#'
#' @description Create a kernel object, to use as variable in a model formula.
#'
#' @param x a formula, a vector or a matrix of variables grouped in the same kernel. It could also be a symetric matrix representing the Gram matrix, associated to a kernel function, already computed by the user.
#' @param kernel.function type of kernel. Possible values are \code{"gaussian"},  \code{"linear"}, \code{"polynomial"}, \code{"sigmoid"}, \code{"inverse.quadratic"} or \code{"equality"}. See details below. If \code{x} is a Gram matrix, associated to a kernel function, already computed by the user, \code{kernel.function} should be equal to \code{"gram.matrix"}.
#' @param scale boolean indicating if variables should be scaled before computing the kernel.
#' @param rho,gamma,d kernel function hyperparameters. See details below.
#'
#' @details To use inside kspm() function. Given two \eqn{p-}dimensional vectors \eqn{x} and \eqn{y},
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
#' @importFrom stats as.formula model.frame na.pass
#'
#'
#' @return A Kernel object including all parameters needed in computation of the model
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @references Liu, D., Lin, X., and Ghosh, D. (2007). Semiparametric regression of multidimensional genetic pathway data: least squares kernel machines and linear mixed models. Biometrics, 63(4), 1079:1088.
#'
#'
#'
#' @export



Kernel <- function(x, kernel.function, scale = TRUE, rho = NULL, gamma = NULL, d = NULL)
{

  ################################################################################
  # CALL
  ################################################################################

  call <- match.call(expand.dots = FALSE)
  call.names <- paste0(call)
  call.order <- match(c("x", "kernel.function", "scale", "rho", "gamma", "d"), names(call), 0L)

  ################################################################################
  # SOME CHECKS
  ################################################################################

  if (is.null(kernel.function)) {
    stop("kernel.function is missing")
  }

  # parameters needed for gaussian kernel
  if (kernel.function == "gaussian" & !is.null(rho)) {
    if (rho <= 0) {
      stop("rho should be positive")
    }
  }
  if (kernel.function == "gaussian" & !is.null(gamma)) {
    gamma <- NULL
    warning("gamma = ", call.names[call.order[5]], " is not used with gaussian kernel")
  }
  if (kernel.function == "gaussian" & !is.null(d)) {
    d <- NULL
    warning("d = ", call.names[call.order[6]], " is not used with gaussian kernel")
  }
  # parameters needed for linear kernel
  if (kernel.function == "linear" & !is.null(rho)) {
    rho <- NULL
    warning("rho = ", call.names[call.order[4]], " is not used with linear kernel")
  }
  if (kernel.function == "linear" & !is.null(gamma)) {
    gamma <- NULL
    warning("gamma = ", call.names[call.order[5]], " is not used with linear kernel")
  }
  if (kernel.function == "linear" & !is.null(d)) {
    d <- NULL
    warning("d = ", call.names[call.order[6]], " is not used with linear kernel")
  }
  # parameters needed for polynomial kernel
  if (kernel.function == "polynomial" & !is.null(rho)) {
    if (rho <= 0) {
      stop("rho should be positive")
    }
  }
  if (kernel.function == "polynomial" & !is.null(d)) {
    if (kernel.function == "polynomial" & (d <= 0 | !check.integer(d))) {
      stop("d should be an integer > 0")
    }
  }
  if (kernel.function == "polynomial" & is.null(d)) {
    stop("d is missing")
  }
  # parameters needed for sigmoid kernel
  if (kernel.function == "sigmoid" & !is.null(d)) {
    d <- NULL
    warning("d = ", call.names[call.order[6]], " is not used with sigmoid kernel")
  }
  # parameters needed for inverse quadratic kernel
  if (kernel.function == "inverse.quadratic" & !is.null(rho)) {
    rho <- NULL
    warning("rho = ", call.names[call.order[4]], " is not used with inverse quadratic kernel")
  }
  if (kernel.function == "inverse.quadratic" & !is.null(gamma)) {
    if (gamma < 0) {
      stop("gamma should be positive")
    }
  }
  if (kernel.function == "inverse.quadratic" & !is.null(d)) {
    d <- NULL
    warning("d = ", call.names[call.order[6]], " is not used with inverse quadratic kernel")
  }
  # parameters needed for equality kernel
  if (kernel.function == "equality" & !is.null(rho)) {
    rho <- NULL
    warning("rho = ", call.names[call.order[4]], " is not used with equality kernel")
  }
  if (kernel.function == "equality" & !is.null(gamma)) {
    gamma <- NULL
    warning("gamma = ", call.names[call.order[5]], " is not used with equality kernel")
  }
  if (kernel.function == "equality" & !is.null(d)) {
    d <- NULL
    warning("d = ", call.names[call.order[6]], " is not used with equality kernel")
  }
  # parameters needed for user kernel
  if (kernel.function == "gram.matrix" & !is.null(rho)) {
    rho <- NULL
    warning("rho = ", call.names[call.order[4]], " is not used with gram.matrix kernel")
  }
  if (kernel.function == "gram.matrix" & !is.null(gamma)) {
    gamma <- NULL
    warning("gamma = ", call.names[call.order[5]], " is not used with gram.matrix kernel")
  }
  if (kernel.function == "gram.matrix" & !is.null(d)) {
    d <- NULL
    warning("d = ", call.names[call.order[6]], " is not used with gram.matrix kernel")
  }

  ################################################################################
  # COMPUTE KERNEL
  ################################################################################

  # if x is a formula
  if (inherits(x, "formula")) {
    kernel.formula <- x
    Z <- model.frame(x, na.action = na.pass)
    K <- NULL
    type.object <- "formula"
  }

  # if x is a vector
  if (is.vector(x)) {
    kernel.formula <- as.formula(paste0("~", call.names[call.order[1]]))
    Z <- matrix(x, nrow = length(x))
    colnames(Z) <- call.names[call.order[1]]
    K <- NULL
    type.object <- "vector"
  }

  # if x is a matrix
  if (is.matrix(x)) {
    if (kernel.function != "gram.matrix") {
      K <- NULL
      # case where x is a symmetric matrix (should be gram.matrix)
      if (isSymmetric(x)) {
        warning("kernel is a symmetric matrix and kernel.function differs from gram.matrix")
      }
      if (!is.null(colnames(x))) {
        Z <- x
        kernel.formula <- as.formula(paste0("~", paste(colnames(Z), collapse = " + ")))
      }
      if (is.null(colnames(x))) {
        Z <- x
        colnames(Z) <- paste0(call.names[call.order[1]], 1:(dim(Z)[2]))
        kernel.formula <- as.formula(paste0("~", paste(colnames(Z), collapse = " + ")))
      }
      type.object <- "design matrix"
    } else {
      if (!isSymmetric(x)) {
        stop("kernel should be a symmetric matrix when kernel.function = gram.matrix")
      }
      K <- x
      Z <- matrix(NA, nrow = dim(K)[1], ncol = 0)
      kernel.formula <- NA
      type.object <- "Gram matrix"
    }
  }

  ################################################################################
  # RETURN OUTPUT
  ################################################################################

  return(list(K = K, Z = Z, kernel.function = kernel.function, kernel.formula = kernel.formula, rho = rho, gamma = gamma, d = d, kernel.scale = scale, type.object = type.object))
}

