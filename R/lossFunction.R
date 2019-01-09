#' @title Computation of the leave one out error (LOOE) in kernel semi parametric model
#'
#' @description internal function to optimize model for estimating hyperparameters based on LOOE
#'
#' @param param. initial parameter values.
#' @param Y. response matrix.
#' @param X. X matrix (linear part).
#' @param kernelList. list of kernels (kernel part).
#' @param n. nb of samples.
#' @param not.missing. nb of non missing samples.
#' @param compute.kernel. boolean. If TRUE, the kernel matrix is computed at each iteration. Should be TRUE when hyperparameters of kernel functions should be estimated by the model.
#' @param print.lambda. boolean. If TRUE, values of tunning parameters (lambda) are printed at each iteration.
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @export




lossFunction.looe <- function(param. = NULL, Y. = NULL, X. = NULL, kernelList. = NULL, n. = NULL, not.missing. = NULL, compute.kernel. = NULL, print.lambda. = FALSE){
  if (print.lambda.) {
    print(param.)
  }
  return(get.parameters( X = X., Y = Y., kernelList = kernelList., free.parameters = param., n = n., not.missing = not.missing., compute.kernel = compute.kernel.)[[1]])
}



