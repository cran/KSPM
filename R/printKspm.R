#' @title Print results from a Kernel Semi parametric Model Fit
#'
#' @description print method for class "kspm".
#'
#'
#'
#' @param x an object used to select a method. Usually, a result of a call to \code{kspm} or a result from \code{summary.kspm}.
#' @param ... additional optional argument (currently unused).
#'
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kspm} for fitting model, \link{summary.kspm}
#'
#'
#' @rdname print.kspm
#' @export print.kspm
#' @export

print.kspm <- function(x, ...)
{
  #-----------------------------------------------------------------------------
  # CALL
  #-----------------------------------------------------------------------------
  cat("\nCall:\n")
  base::print(x$call)

  #-----------------------------------------------------------------------------
  # COEFFICIENTS
  #-----------------------------------------------------------------------------
  cat("\nLinear coefficients:\n")
  base::print(round(x$linear.coefficients[,1], 5))
}


#' @rdname print.kspm
#' @export print.summary.kspm
#' @export

print.summary.kspm <- function(x, ...) {

  #-----------------------------------------------------------------------------
  # CALL
  #-----------------------------------------------------------------------------
  cat("\nCall:\n")
  base::print(x$call)

  #-----------------------------------------------------------------------------
  # SAMPLE SIZE
  #-----------------------------------------------------------------------------
  cat("\nSample size:")
  cat(paste0("\nn = ", x$sample.size$inc))
  # Indicate number of missing data that were removed
  if (x$sample.size$inc != x$sample.size$all) {
    cat(paste0("\n(", x$sample.size$all - x$sample.size$inc, " observation deleted due to missingness)"))
  }

  #-----------------------------------------------------------------------------
  # SUMMARY OF RESIDUALS
  #-----------------------------------------------------------------------------
  cat("\n\nResiduals: \n")
  summaryResiduals <- round(base::summary(x$residuals)[c(1:3,5,6)], 4)
  names(summaryResiduals) <- c("Min", "Q1", "Median", "Q3", "Max")
  base::print(summaryResiduals)


  #-----------------------------------------------------------------------------
  # LINEAR PART
  #-----------------------------------------------------------------------------
  if (!is.null(x$coefficients)) {
    cat("\nCoefficients (linear part): \n")
    base::print(x$coefficients)
  } else {
    cat("\nNo coefficient in linear part \n")
  }

  #-----------------------------------------------------------------------------
  # KERNEL PART
  #-----------------------------------------------------------------------------
  if (!is.null(x$score.test)) {
    cat("\nScore test for non-parametric kernel: \n")
    base::print(noquote(x$score.test))
  }
  if (!is.null(x$global.p.value)) {
    cat(paste0("\nGlobal test: p-value = ", round(x$global.p.value, 4),"\n"))
  }


  #-----------------------------------------------------------------------------
  # EDF AND R2/R2adj
  #-----------------------------------------------------------------------------
  cat(paste0("\nResidual standard error: ", round(x$sigma, 2), " on ", round(x$edf, 2), " effective degrees of freedom"))
  cat(paste0("\nMultiple R-squared: ", round(x$R2, 4), ", Adjusted R-squared: ", round(x$R2adj, 4)))
  cat("\n ")

}



