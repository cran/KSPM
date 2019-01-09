#' @title Extract AIC from a Kernel Semi Parametric Model
#'
#' @description Computes the Akaike Information Criterion (AIC) for a kspm fit.
#'
#' @param fit fitted model, usually the result of \link{kspm}.
#' @param scale option not available for kspm fit.
#' @param k numeric specifying the 'weight' of the effective degrees of freedom (edf) part in the AIC formula. See details.
#' @param correction boolean indicating if the corrected AIC should be computed instead of standard AIC, may be \code{TRUE} only for \code{k=2}. See details.
#' @param ... additional optional argument (currently unused).
#'
#' @details The criterion used is \eqn{AIC = n log(RSS) + k (n-edf)}{AIC = n log(RSS) + k (n-edf)} where \eqn{RSS}{RSS} is the residual sum of squares and \eqn{edf}{edf} is the effective degree of freedom of the model. \code{k = 2} corresponds to the traditional AIC, using \code{k = log(n)} provides Bayesian Information Criterion (BIC) instead. For \code{k=2}, the corrected Akaike's Information Criterion (AICc) is obtained by \eqn{AICc = AIC + \frac{2 (n-edf) (n-edf+1)}{(edf-1)}}{AICc = AIC + 2*(n-edf)*(n-edf+1) / (edf-1)}.
#'
#' @return \code{extractAIC.kspm} returns a numeric value corresponding to AIC. Of note, the AIC obtained here differs from a constant to the AIC obtained with \code{extractAIC} applied to a \link{lm} object. If one wants to compare a \code{kspm} model with a \code{lm} model, it is preferrable to compute again the \code{lm} model using \link{kspm} function by specifying \code{kernel = NULL} and apply \code{extractAIC} method on this model.
#'
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @references Liu, D., Lin, X., and Ghosh, D. (2007). Semiparametric regression of multidimensional genetic pathway data: least squares kernel machines and linear mixed models. Biometrics, 63(4), 1079:1088.
#'
#' @seealso \link{stepKSPM} for variable selection procedure based on AIC.
#'
#' @examples
#' x <- 1:15
#' y <- 3*x + rnorm(15, 0, 2)
#' fit <- kspm(y, kernel = ~ Kernel(x, kernel.function = "linear"))
#' extractAIC(fit)
#'
#' @importFrom stats extractAIC
#'
#' @rdname extractAIC.kspm
#' @export extractAIC.kspm
#' @export






extractAIC.kspm = function(fit, scale = NULL, k = 2, correction = FALSE, ...){

  ################################################################################
  # SOME CHECKS
  ################################################################################

  # Warning if correction is asking for k different from 2
  if (correction & k != 2) {
    warning("correction is not applied (valid only for k = 2), see ?extractAIC.kspm for more details")
  }

  ################################################################################
  # Computation
  ################################################################################

  # Compute the Residual Sum of Squares
  RSS <- (t(fit$residuals) %*% (fit$residuals))[1, 1]
  # Compute AIC
  out <- fit$n * log(RSS) + k * (fit$n - fit$edf)

  # Compute AICc
  if (correction & k == 2) {
    out.correction = out + 2 * (fit$n - fit$edf) * (fit$n - fit$edf + 1) / (fit$edf - 1)
    return(out.correction)
  } else {
    return(out)
  }

}



