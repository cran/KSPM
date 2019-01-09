#' @title Plot derivatives of a kspm object
#'
#' @description Plot of derivatives for kernel part of a kspm model.
#'
#'
#' @param x an object of class "derivatives", usually, a result of a call to \code{derivatives}.
#' @param subset if a subset of the plots is required, specify the names of the variable for which plot of derivatives is required.
#' @param xlab x label
#' @param ylab y label
#' @param ... further arguments passed to or from other methods.
#'
#' @details X axis represents the raw data used as input in kernel part of the model. Y axis represents the pointwise derivative values i.e. the derivatives of fitted value according to the variable of interest.
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @seealso \link{derivatives}
#'
#' @examples
#' x <- 1:15
#' z1 <- runif(15, 1, 6)
#' z2 <- rnorm(15, 1, 2)
#' y <- 3*x + (z1 + z2)^2 + rnorm(15, 0, 2)
#' fit <- kspm(y, linear = ~ x, kernel = ~ Kernel(~ z1 + z2,
#' kernel.function = "polynomial", d= 2, rho = 1, gamma = 0))
#' plot(derivatives(fit))
#'
#' @references
#'
#' Kim, Choongrak, Byeong U. Park, and Woochul Kim. "Influence diagnostics in semiparametric regression models." Statistics and probability letters 60.1 (2002): 49:58.
#'
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics abline legend lines plot
#'
#' @export
#'


plot.derivatives <- function(x, subset = NULL, xlab = NULL, ylab = NULL, ...)
{

  ################################################################################
  # CHECKS
  ################################################################################

  # object should be of class 'derivatives'
  if (!inherits(x, "derivatives")) {
    warning("calling plot.derivatives(<fake-derivatives-object>) ...")
  }


  ################################################################################
  # Z SUBSET OF VARIABLES FOR WHICH GRAPHICS SHOULD BE DISPLAYED
  ################################################################################

  # if NULL (default), a graphic is displayed for each variable
  if (is.null(subset)) {
    subset <- rownames(x$modelmat)
  } else {
    if (!is.vector(subset) | !is.character(subset)) {
      stop("invalid argument for subset, see ?plot.derivatives for more details")
    }
    if (sum(!(subset %in% colnames(x$rawmat))) > 0) {
      stop("subset contains variable(s) not used in the model")
    }
  }


  ################################################################################
  # PLOT
  ################################################################################

  # one graphic by variable for which a graphic was asked
  for (variable.name in subset) {

    #-------------------------------------------------------------------------------
    # X-axis
    #-------------------------------------------------------------------------------
    x1 <- x$rawmat[, variable.name]

    #-------------------------------------------------------------------------------
    # Y-axis
    #-------------------------------------------------------------------------------

    # Number of times the variable appears in a kernel
    yn <- sum(x$modelmat[variable.name, ], na.rm = TRUE)
    # Indicates if the variable appears in a kernel that does not represent an interaction
    yn2 <- x$modelmat[variable.name, names(x$derivmat)]
    names(yn2) <- names(x$derivmat)
    # Number of times the variable appears in a kernel that does not represent an interaction
    yn3 <- sum(x$modelmat[variable.name, names(x$derivmat)], na.rm = TRUE)
    kernelAppearence <- na.omit(names(yn2[yn2 == 1])) # kernel in which the variable appear

    #-------------------------------------------------------------------------------
    # Plot if variable appears only in one kernel (not interaction)
    #-------------------------------------------------------------------------------

    if (yn3 == 1) {
      # set the value of xlab and ylab
      if (is.null(xlab)) {
        xlab.current <- paste0(variable.name, " (", kernelAppearence, ")")
      } else {
        xlab.current <- xlab
      }

      if (is.null(ylab)) {
        ylab.current <- paste("Derivatives dh/d", variable.name, " (", kernelAppearence, ")")
      } else {
        ylab.current <- ylab
      }
      plot(x1, x$derivmat[[kernelAppearence]][, variable.name], xlab = xlab.current, ylab = ylab.current, ...)
      abline(h = 0, col = "gray")
      # Warning if the variable also appear in interaction
      if (yn3 != yn) {
        legend("bottomleft", legend = c("warning: main effects only", "interactions are not displayed"), bty = "n")
      }
    } else {

      #-------------------------------------------------------------------------------
      # Plot if variable appears in several kernels (not interaction)
      #-------------------------------------------------------------------------------

      # list of values of derivatives for each kernel
      y <- list()
      for (i in 1:length(kernelAppearence)) {
        y[[i]] <- x$derivmat[[kernelAppearence[i]]][, variable.name]

      }
      yT <- Reduce("+", y)
      ymin <- Reduce("min", y)
      ymax <- Reduce("max", y)

      # set the value of xlab and ylab
      if (is.null(xlab)) {
        xlab.current <- variable.name
      } else {
        xlab.current <- xlab
      }
      if (is.null(ylab)) {
        ylab.current <- paste("Derivatives dh/d", variable.name)
      } else {
        ylab.current <- ylab
      }
      plot(x1, yT, ylim = c(ymin, ymax), xlab = xlab.current, ylab = ylab.current, ...)
      # add sub value for each kernel
      for (i in 1:length(kernelAppearence)) {
        lines(x1, y[[i]], type = "p", col = i)
      }
      legend("bottomright", legend = c(kernelAppearence, "All"), col = c(1:length(kernelAppearence), 1), bty = "n", pch = c(rep(1, length(kernelAppearence)), 16))
      abline(h = 0, col = "gray")
      if (yn3 != yn) {
        legend("bottomleft", legend = c("warning: main effects only", "interactions are not displayed"), bty = "n")
      }
    }
    # ask for new graphic
    devAskNewPage(ask = TRUE)

  }

  # stop asking for new graphic
  devAskNewPage(ask = FALSE)

}


