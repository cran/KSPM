#' @title Giving information about Kernel Semi parametric Model Fits
#'
#' @description gives information about Kernel Semi parametric Model Fits
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param print logical, if \code{TRUE}, table of information are printed.
#'
#'
#' @return \code{info.kspm} returns a table of information whose each row corresponds to a kernel included in the model and columns are:
#'   \item{type}{type of object used to define the kernel}
#'   \item{dim}{dimension of data used in the model}
#'   \item{type.predict}{type of object the user should provide in \link{predict.kspm} function}
#'   \item{dim.predict}{dimension of object the user should provide in \link{predict.kspm} function}
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#'
#' @seealso \link{kspm}, \link{predict.kspm}
#'
#' @export
#'


info.kspm <- function(object, print = TRUE)
{


  ################################################################################
  # SOME CHECKS
  ################################################################################

  # object should be of class kspm
  if (!inherits(object, "kspm")) {
    stop("object should be of class 'kspm'")
  }

  ################################################################################
  # COMPUTE INFORMATION
  ################################################################################

  # empty data.frame
  out <- data.frame(type = NA, dim = NA, type.predict = NA, dim.predict = NA)

  for (k in 1:length(names(object$kernel.info))) {
    out[k, "type"] <- object$kernel.info[[k]]$type.object
    if (object$kernel.info[[k]]$type.object == "formula") {
      out[k, "dim"] <- paste0("n=", object$n, " x ", "q=", dim(object$kernel.info[[k]]$Z)[2])
      out[k, "type.predict"] <- "design matrix or data.frame with column names"
      out[k, "dim.predict"] <- paste("nn", "x", dim(object$kernel.info[[k]]$Z)[2])
    }
    if (object$kernel.info[[k]]$type.object == "vector") {
      out[k, "dim"] <- paste(object$n)
      out[k, "type.predict"] <- "vector or design matrix with 1 column"
      out[k, "dim.predict"] <- "nn or nn x 1"
    }
    if (object$kernel.info[[k]]$type.object == "design matrix") {
      out[k, "dim"] <- paste0("n=", object$n, " x ", "q=", dim(object$kernel.info[[k]]$Z)[2])
      out[k, "type.predict"] <- "design matrix (same order of columns)"
      out[k, "dim.predict"] <- paste("nn", "x", dim(object$kernel.info[[k]]$Z)[2])
    }
    if (object$kernel.info[[k]]$type.object == "Gram matrix") {
      out[k, "dim"] <- paste0("n=", object$n, " x ", "n=", object$n)
      out[k, "type.predict"] <- "Gram matrix for similarity of nn samples with n samples"
      out[k, "dim.predict"] <- paste("nn", "x", object$n)
    }
  }

  # print results
  row.names(out) <- names(object$kernel.info)
  print(out)
  # return results
  return(invisible(out))
}
