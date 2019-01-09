#' @title Control various aspects of the optimisation problem
#'
#' @description Allow the user to set some characteristics of the optimisation algorithm
#'
#' @param interval.upper integer or vetor of initial maximum value(s) allowed for parameter(s)
#' @param interval.lower integer or vetor of initial maximum value(s) allowed for parameter(s)
#' @param trace boolean. If TRUE parameters value at each iteration are displayed.
#' @param optimize.tol if \link[stats]{optimize} function is used. See \link[stats]{optimize}
#' @param NP if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#' @param itermax if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#' @param CR if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#' @param F if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#' @param initialpop if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#' @param storepopfrom if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#' @param storepopfreq if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#' @param p if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#' @param c if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#' @param reltol if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#' @param steptol if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#' @param parallel if \link[DEoptim]{DEoptim} function is used. See \link[DEoptim]{DEoptim.control}
#'
#' @details When only one hyperparameter should be estimated, the optimisation problem calls the \link[stats]{optimize} function from \code{stats} basic package. Otherwise, it calls the \link[DEoptim]{DEoptim} function from the package \code{DEoptim}. In both case, the parameters are choosen among the initial interval defined by \code{interval.lower} and \code{interval.upper}.
#'
#' @return \code{search.parameters} is an iterative algorithm estimating model parameters and returns the following components:
#'  \item{lambda}{tuning parameters for penalization.}
#'  \item{beta}{vector of coefficients associated with linear part of the model, the size being the number of variable in linear part (including an intercept term).}
#'  \item{alpha}{vector of coefficients associated with kernel part of the model, the size being the sample size.}
#'  \item{Ginv}{a matrix used in several calculations. \eqn{Ginv = (\lambda I + K)^{-1}}{Ginv = (lambda I + K)^(-1)}.}
#'
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso link get.parameters for computation of parameters at each iteration
#'
#'
#' @export



kspmControl <- function(
  # KSPM
  interval.upper = NA, interval.lower = NA, trace = FALSE,
  # Optimize
  optimize.tol = .Machine$double.eps^0.25,
  # DEoptim
  NP = NA, itermax = 500, CR = 0.5, F = 0.8, initialpop = NULL,
  storepopfrom = itermax + 1, storepopfreq = 1, p = 0.2, c = 0, reltol = sqrt(.Machine$double.eps),
  steptol = itermax,
  # parallel in DEoptim
  parallel = FALSE)
{

  return(list(interval.upper = interval.upper, interval.lower = interval.lower,

              # Optimize
              optimize.tol = optimize.tol,

              # DEoptim
              NP = NP, itermax = itermax, CR = CR, F = F, trace = trace, initialpop = initialpop,
              storepopfrom = storepopfrom, storepopfreq = storepopfreq, p = p, c = c, reltol = reltol,
              steptol = steptol,

              # parallele in DEoptim
              parallel = parallel))
}

