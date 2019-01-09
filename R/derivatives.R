#' @title Computing kernel function derivatives
#'
#' @description \code{derivatives} is a function for "kspm" object computing pointwise partial derivatives of \eqn{h(Z)} accroding to each \eqn{Z} variable.
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#'
#'
#' @details derivatives are not computed for interactions. If a variable is included in several kernels, the user may obtain the corresponding pointwise derivatives by summing the pointwise derivatives associated with each kernel.
#'
#'
#' @return an object of class 'derivatives'
#' \item{derivmat}{a list of \eqn{n \times d}{n x d} matrix (one for each kernel) where \eqn{n}{n} is the number of subjects and \eqn{d}{d} the number of variables included in the kernel}
#' \item{rawmat}{a \eqn{n \times q}{n x q} matrix with all variables included in the kernel part of the model \eqn{q}{q} the number of variables included in the whole kernel part}
#' \item{scalemat}{scaled version of rawmat}
#' \item{modelmat}{matrix of correspondance between variable and kernels}
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @seealso \link{plot.derivatives}
#'
#' @references
#'
#' Kim, Choongrak, Byeong U. Park, and Woochul Kim. "Influence diagnostics in semiparametric regression models." Statistics and probability letters 60.1 (2002): 49:58.
#'
#'
#' @export
#'



derivatives <- function(object)
{

  ################################################################################
  # SOME CHECKS
  ################################################################################

  # object should be of class kspm
  if (!inherits(object, "kspm")) {
    stop("object should be of class 'kspm'")
  }


  ################################################################################
  # USEFUL PARAMETERS
  ################################################################################

  # sample size
  n <- object$n

  # list of matrix of derivatives (one by kernel)
  derivmatL <- list()
  # list of which variable in which kernel
  listVariableL <- c()

  # design matrix for variables in the kernel
  row <- matrix(NA, nrow = n, ncol = 0)
  rownames(row) <- rownames(object$K[[1]])

  # design matrix for variables in the kernel, after scaling
  row.scale <- matrix(NA, nrow = n, ncol = 0)
  rownames(row.scale) <- rownames(object$K[[1]])


  ################################################################################
  # MATRIX OF DERIVATIVES
  ################################################################################

  # a matrix is computed for each kernel (interaction excluded)
  for (l in 1:length(object$kernel.info)) {
    if (object$kernel.info[[l]]$kernel.function == "gram.matrix") {
      warning(paste0("Derivatives cannot be computed for ", names(object$kernel.info)[l], " because it was defined using Gram matrix"))
      next
    }
    # kernel name
    nameKernel <- names(object$kernel.info)[l]

    # kernel type
    whichkernel <- object$kernel.info[[l]]$kernel.function

    # Z matrix
    if (object$kernel.info[[l]]$kernel.scale) {
      Z2 <- object$kernel.info[[l]]$Z.scale
    } else {
      Z2 <- object$kernel.info[[l]]$Z
    }

    # Z dimension
    q <- dim(Z2)[2]

    # create empty matrix to record derivatives: 1 value by sample (n rows) and by variable in Kernel function (q named columns)
    derivmat <- matrix(NA, n, q)
    colnames(derivmat) <- colnames(Z2)


    #-------------------------------------------------------------------------------
    # Compute derivative in case of Gaussian kernel
    #-------------------------------------------------------------------------------
    if (whichkernel == "gaussian") {
      rows <- cbind(rep(1:n, each = n), 1:n)
      distances <- Z2[rows[, 1], ] - Z2[rows[, 2], ]
      for (k in 1:q) {
        if (q == 1) {
          distk <- matrix(distances, n, n, byrow = TRUE)
        } else{
          distk <- matrix(distances[, k], n, n, byrow = TRUE)
        }
        derivmat[, k] <- (-2/object$kernel.info[[l]]$rho) * (distk * object$K[nameKernel][[1]]) %*% object$kernel.coefficients[, nameKernel] # pointwise derivatives
      }
    }

    #-------------------------------------------------------------------------------
    # Compute derivative in case of Linear kernel
    #-------------------------------------------------------------------------------
    if (whichkernel == "linear") {
      for (k in 1:q) {
        derivmat[, k] <- Z2[, k] %*% object$kernel.coefficients[, nameKernel] # pointwise derivatives
      }
    }

    #-------------------------------------------------------------------------------
    # Compute derivative in case of Polynomial kernel
    #-------------------------------------------------------------------------------
    if (whichkernel == "polynomial") {
      K_d_1 = (object$kernel.info[[l]]$rho * Z2 %*% t(Z2) + object$kernel.info[[l]]$gamma) ^ (object$kernel.info[[l]]$d - 1)
      for (k in 1:q) {
        derivmat[, k] <- object$kernel.info[[l]]$d * (K_d_1 %*% object$kernel.coefficients[, nameKernel] %*% Z2[, k])[ ,k]  # pointwise derivatives
      }
    }

    #-------------------------------------------------------------------------------
    # Compute derivative in case of sigmoidal kernel
    #-------------------------------------------------------------------------------
    if (whichkernel == "sigmoid") {
      for (k in 1:q) {
        derivmat[, k] <- object$kernel.info[[l]]$rho * ((1 - (object$K[, nameKernel])^2) %*% object$kernel.coefficients[, nameKernel] %*% Z2[, k] )[ ,k]  # pointwise derivatives
      }
    }

    #-------------------------------------------------------------------------------
    # Compute derivative in case of inverse quadratic kernel
    #-------------------------------------------------------------------------------
    if (whichkernel == "inverse.quadratic") {
      rows <- cbind(rep(1:n, each = n), 1:n)
      distances <- Z2[rows[, 1], ] - Z2[rows[, 2], ]
      for (k in 1:q) {
        if (q == 1) {
          distk <- matrix(distances, n, n, byrow = TRUE)
        } else{
          distk <- matrix(distances[, k], n, n, byrow = TRUE)
        }
        derivmat[, k] <- -(distk * (object$K[, nameKernel])^3) %*% object$kernel.coefficients[, nameKernel] # pointwise derivatives
      }
    }

    #-------------------------------------------------------------------------------
    # case of kernel defined by user.matrix or equality kernel
    #-------------------------------------------------------------------------------
    if (whichkernel %in% c("user.matrix", "equality")) {
      derivmat <- NULL
    }

    # save the matrix in the list
    derivmatL[[nameKernel]] <- derivmat
    # save the names of the variables as associated with this kernel
    listVariableL <- c(listVariableL, colnames(derivmat))
  }


  ################################################################################
  # WHICH VARIABLE IN WHICH KERNEL ?
  ################################################################################

  # built the matrix in which information will take place
  listVariable <- unique(listVariableL)
  design <- matrix(NA, nrow = length(listVariable), ncol = length(object$lambda))
  # rows are all possible variable, columns are all possible kernel part
  colnames(design) <- names(object$lambda)
  row.names(design) <- c(listVariable)

  # for each kernel part (interaction included)
  for (l in 1:length(object$lambda)) {

    # name of main kernels composing this kernel part
    nom <- strsplit(names(object$lambda)[l], ":")[[1]]

    # for each main kernel
    for (i in 1:length(nom)) {

      # variables in this main kernel
      colnom <- colnames(object$kernel.info[[nom[i]]]$Z)
      # fill matrix of information about variables that are in this kernel
      design[colnom, names(object$lambda)[l]] <- 1

      # for each variables in this main kernel, we add the variable in row and row.scale if it is not already included
      if (length(colnom) != 0) {
        for (j in 1:length(colnom)) {
          if (!(colnom[j] %in% colnames(row))) {
            row2 <- cbind(row, object$kernel.info[[nom[i]]]$Z[, colnom[j]])
            colnames(row2) <- c(colnames(row), colnom[j])
            row <- row2
            if (!is.na(object$kernel.info[[nom[i]]]$Z.scale[1])) {
              row.scale2 <- cbind(row.scale, object$kernel.info[[nom[i]]]$Z.scale[, colnom[j]])
            } else {
              row.scale2 <- cbind(row.scale, matrix(rep(NA, n), ncol = 1))
            }
            colnames(row.scale2) <- c(colnames(row.scale), colnom[j])
            row.scale <- row.scale2
          }
        }
      }
    }
  }


  ################################################################################
  # RETURN OUTPUT
  ################################################################################

  # Return an object of class derivatives
  out <- list(derivmat = derivmatL, rawmat = row, scalemat = row.scale, modelmat = design)
  class(out) <- c("derivatives")
  return(out)

}


