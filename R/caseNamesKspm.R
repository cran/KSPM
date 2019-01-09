#' @title Case names of fitted models
#'
#' @description Simple utility returning names of cases involved in a kernel semi parametric model.
#'
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param ... additional optional argument (currently unused).
#'
#'
#'
#' @return a character vector.
#'
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kspm} for fitting model, \link{nobs.kspm}, \link{variable.names.kspm}.
#'
#' @importFrom stats case.names
#'
#' @rdname case.names.kspm
#' @export case.names.kspm
#' @export


case.names.kspm <- function(object, ...)
{
  return(names(object$fitted.values))
}



