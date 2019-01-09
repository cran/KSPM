#' @title Optimisation to cumpute hyperparameter in Kernel Semi Parametric model
#'
#' @description internal function to optimize model for estimating hyperparameters
#'
#' @param Y response matrix
#' @param X X matrix
#' @param kernelList of kernels
#' @param n nb of samples
#' @param not.missing nb of non missing samples
#' @param compute.kernel boolean kernel computation
#' @param controlKspm control parameters
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @importFrom stats optimize
#' @importFrom DEoptim DEoptim DEoptim.control
#'
#' @export






# search.parameters <- function(Y = NULL, X = NULL, kernelList = NULL,  n = NULL, not.missing = NULL, compute.kernel = NULL, controlKspm = NULL)
# {
#
#   ################################################################################
#   # DEFINE CONTROL PARAMETERS 'interval.upper' AND 'interval.lower'
#   ################################################################################
#
#   # set an initial interval.upper value if it has not been provided by the user
#   if (is.na(controlKspm$interval.upper)) {
#     if (compute.kernel) {
#       controlKspm$interval.upper <- rep(n, length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel))
#     } else {
#       controlKspm$interval.upper <- rep(n, length(kernelList))
#     }
#   } else {# if interval.upper is provided, it should respect constraints
#     if (compute.kernel) {
#       # check on length
#       if (length(controlKspm$interval.upper) != length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel)) {
#         stop(paste0("'interval.upper' should be of length ", length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel)))
#       }
#       # check on values
#       test.constraint <- sum(controlKspm$interval.upper <= rep(Inf, length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel)))
#     } else {
#       # check on length
#       if (length(controlKspm$interval.upper) != length(kernelList$inclKernel)) {
#         stop(paste0("'interval.upper' should be of length ", length(kernelList$inclKernel)))
#       }
#       # check on values
#       test.constraint <- sum(controlKspm$interval.upper <= rep(Inf, length(kernelList$inclKernel)))
#     }
#     if (test.constraint > 0) {
#       stop("initial value for 'interval.upper' does not respect constraints of hyperparameters")
#     }
#   }
#
#   # set an initial interval.lower value if it has not been provided by the user
#   if (is.na(controlKspm$interval.lower)) {
#     # the interval.lower values depends on tuning parameter constraints
#     # for lambda parameters the constraint is that they should be positive
#     if (compute.kernel) {
#       controlKspm$interval.lower <- c(kernelList$free.parameters.contraint, rep(0, length(kernelList$inclKernel)))
#     } else {
#       controlKspm$interval.lower <- rep(0, length(kernelList))
#     }
#   } else {# if interval.lower is provided, it should respect constraints
#     if (compute.kernel) {
#       # check on length
#       if (length(controlKspm$interval.lower) != length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel)) {
#         stop(paste0("'interval.lower' should be of length ", length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel)))
#       }
#       # check on values
#       test.constraint <- sum(controlKspm$interval.lower >= c(kernelList$free.parameters.contraint, rep(0, length(kernelList$inclKernel))))
#     } else {
#       # check on length
#       if (length(controlKspm$interval.lower) != length(kernelList$inclKernel)) {
#         stop(paste0("'interval.lower' should be of length ", length(kernelList$inclKernel)))
#       }
#       # check on values
#       test.constraint <- sum(controlKspm$interval.lower >= rep(0, length(kernelList$inclKernel)))
#     }
#     if (test.constraint > 0) {
#       stop("initial value for 'interval.lower' does not respect constraints of hyperparameters")
#     }
#   }
#
#
#   ################################################################################
#   # CASE OF ONE KERNEL WITHOUT TUNING PARAMETER
#   ################################################################################
#
#   if (length(controlKspm$interval.upper) == 1) {
#     lambdaOpt <- optimize(f = lossFunction.looe, interval = c(controlKspm$interval.lower, controlKspm$interval.upper), Y. = Y, X. = X, kernelList. = kernelList, n. = n, not.missing. = not.missing, compute.kernel. = compute.kernel, print.lambda. = controlKspm$trace, tol = controlKspm$optimize.tol)$minimum
#     # final parameters
#     out <- get.parameters(Y = Y, X = X, kernelList = kernelList, free.parameters = lambdaOpt, n = n, not.missing = not.missing, compute.kernel = compute.kernel)
#   }
#
#
#   ################################################################################
#   # OTHER CASES
#   ################################################################################
#
#   if (length(controlKspm$interval.lower) > 1) {
#     # DEoptim does not accept 0 as border
#     interval.min <- controlKspm$interval.lower + 1e-16
#     interval.max <- controlKspm$interval.upper + 1e-16
#     outDEoptim <- DEoptim(lossFunction.looe, Y. = Y, X. = X, kernelList. = kernelList, n. = n, not.missing. = not.missing, compute.kernel. = compute.kernel, lower = interval.min, upper = interval.max, control = DEoptim.control(trace = controlKspm$trace, itermax = controlKspm$itermax, NP = controlKspm$NP, CR = controlKspm$CR, F = controlKspm$F, initialpop = controlKspm$initialpop, storepopfrom = controlKspm$storepopfrom, storepopfreq = controlKspm$storepopfreq, p = controlKspm$p, c = controlKspm$c, reltol = controlKspm$reltol, steptol = controlKspm$steptol))
#     # best values for hyperparameters
#     lambdaOpt <- outDEoptim$optim$bestmem
#     # give information about convergence
#     # if maximum is reach
#     hyperparameters.name <- c(kernelList$free.parameters, rep("lambda", length(kernelList$inclKernel)))
#     for (i in 1:length(lambdaOpt)) {
#       if (interval.max[i] - lambdaOpt[i] <= 0.001 * interval.max[i]) {
#         warning(paste0("maximum is reached for hyperparameter ", hyperparameters.name[1], ". Initial upper value = ", interval.max[i], ". Estimated value = ", lambdaOpt[i], "."))
#       }
#     }
#     # final parameters
#     out <- get.parameters(Y = Y, X = X, kernelList = kernelList, free.parameters = lambdaOpt, n = n, not.missing = not.missing, compute.kernel = compute.kernel)
#   }
#
#
#   ################################################################################
#   # RETURN PARAMETERS
#   ################################################################################
#
#   return(out)
# }
















search.parameters <- function(Y = NULL, X = NULL, kernelList = NULL,  n = NULL, not.missing = NULL, compute.kernel = NULL, controlKspm = NULL)
{

  ################################################################################
  # DEFINE CONTROL PARAMETERS 'interval.upper' AND 'interval.lower'
  ################################################################################

  # set an initial interval.upper value if it has not been provided by the user
  if (is.na(controlKspm$interval.upper)) {
    if (compute.kernel) {
      controlKspm$interval.upper <- rep(n, length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel))
    } else {

      controlKspm$interval.upper <- rep(n, length(kernelList))

    }
  } else {# if interval.upper is provided, it should respect constraints
    if (compute.kernel) {
      # check on length
      if (length(controlKspm$interval.upper) != length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel)) {
        stop(paste0("'interval.upper' should be of length ", length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel)))
      }
      # check on values
      test.constraint <- sum(controlKspm$interval.upper <= rep(Inf, length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel)))
    } else {

      # check on length
      if (length(controlKspm$interval.upper) != length(kernelList)) {
        stop(paste0("'interval.upper' should be of length ", length(kernelList)))
      }
      # check on values
      test.constraint <- sum(controlKspm$interval.upper <= rep(Inf, length(kernelList)))
    }
    if (test.constraint > 0) {
      stop("initial value for 'interval.upper' does not respect constraints of hyperparameters")
    }
  }

  # set an initial interval.lower value if it has not been provided by the user
  if (is.na(controlKspm$interval.lower)) {
    # the interval.lower values depends on tuning parameter constraints
    # for lambda parameters the constraint is that they should be positive
    if (compute.kernel) {
      controlKspm$interval.lower <- c(kernelList$free.parameters.contraint, rep(0, length(kernelList$inclKernel)))
    } else {
      controlKspm$interval.lower <- rep(0, length(kernelList))
    }
  } else {# if interval.lower is provided, it should respect constraints
    if (compute.kernel) {
      # check on length
      if (length(controlKspm$interval.lower) != length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel)) {
        stop(paste0("'interval.lower' should be of length ", length(kernelList$free.parameters.contraint) + length(kernelList$inclKernel)))
      }
      # check on values
      test.constraint <- sum(controlKspm$interval.lower >= c(kernelList$free.parameters.contraint, rep(0, length(kernelList$inclKernel))))
    } else {
      # check on length
      if (length(controlKspm$interval.lower) != length(kernelList)) {
        stop(paste0("'interval.lower' should be of length ", length(kernelList)))
      }
      # check on values
      test.constraint <- sum(controlKspm$interval.lower >= rep(0, length(kernelList$inclKernel)))
    }
    if (test.constraint > 0) {
      stop("initial value for 'interval.lower' does not respect constraints of hyperparameters")
    }
  }


  ################################################################################
  # CASE OF ONE KERNEL WITHOUT TUNING PARAMETER
  ################################################################################

  if (length(controlKspm$interval.upper) == 1) {

    lambdaOpt <- optimize(f = lossFunction.looe, interval = c(controlKspm$interval.lower, controlKspm$interval.upper), Y. = Y, X. = X, kernelList. = kernelList, n. = n, not.missing. = not.missing, compute.kernel. = compute.kernel, print.lambda. = controlKspm$trace, tol = controlKspm$optimize.tol)$minimum
    # final parameters
    out <- get.parameters(Y = Y, X = X, kernelList = kernelList, free.parameters = lambdaOpt, n = n, not.missing = not.missing, compute.kernel = compute.kernel)
  }


  ################################################################################
  # OTHER CASES
  ################################################################################

  if (length(controlKspm$interval.lower) > 1) {

    # DEoptim does not accept 0 as border
    interval.min <- controlKspm$interval.lower + 1e-16
    interval.max <- controlKspm$interval.upper + 1e-16
    outDEoptim <- DEoptim(lossFunction.looe, Y. = Y, X. = X, kernelList. = kernelList, n. = n, not.missing. = not.missing, compute.kernel. = compute.kernel, lower = interval.min, upper = interval.max, control = DEoptim.control(trace = controlKspm$trace, itermax = controlKspm$itermax, NP = controlKspm$NP, CR = controlKspm$CR, F = controlKspm$F, initialpop = controlKspm$initialpop, storepopfrom = controlKspm$storepopfrom, storepopfreq = controlKspm$storepopfreq, p = controlKspm$p, c = controlKspm$c, reltol = controlKspm$reltol, steptol = controlKspm$steptol, parallelType = ifelse(controlKspm$parallel, 1, 0))) #
    # best values for hyperparameters
    lambdaOpt <- outDEoptim$optim$bestmem
    # give information about convergence
    # if maximum is reach
    hyperparameters.name <- c(kernelList$free.parameters, rep("lambda", length(kernelList$inclKernel)))
    for (i in 1:length(lambdaOpt)) {
      if (interval.max[i] - lambdaOpt[i] <= 0.001 * interval.max[i]) {
        warning(paste0("maximum is reached for hyperparameter ", hyperparameters.name[1], ". Initial upper value = ", interval.max[i], ". Estimated value = ", lambdaOpt[i], "."))
      }
    }
    # final parameters
    out <- get.parameters(Y = Y, X = X, kernelList = kernelList, free.parameters = lambdaOpt, n = n, not.missing = not.missing, compute.kernel = compute.kernel)
  }


  ################################################################################
  # RETURN PARAMETERS
  ################################################################################

  return(out)
}

