#' @title Variable names of fitted models
#'
#' @description Simple utility returning names of variables involved in a kernel semi parametric model.
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param ... additional optional argument (currently unused).
#'
#'
#'
#' @return a list of character vectors. The first element correspond to the names of variables included in the linear part of the model. Then, a vector containing names of variables including in kernel part is provided for each kernel.
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kspm}, \link{summary.kspm}, \link{case.names.kspm}.
#'
#' @importFrom stats variable.names
#'
#'
#' @rdname variable.names.kspm
#' @export variable.names.kspm
#' @export


variable.names.kspm <- function(object, ...)
{
  out <- list(linear = rownames(object$linear.coefficients))
  if (!is.null(object$lambda)) {
    for (i in 1:length(object$lambda)) {
      out[[names(object$lambda)[i]]] <- colnames(object$kernel.info[[names(object$lambda)[i]]]$Z)
    }
  }
  return(out)
}

