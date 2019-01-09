#' @title some internal methods in computation of kernel semi parametric model
#'
#' @description internal methods
#'
#' @param x list of objects
#' @param N numeric value
#' @param object formula provided in the kernel part of \code{kspm} function
#' @param form formula
#' @param sep separator
#' @param ind index value
#' @param nameKernel name of kernel
#' @param not.missing non missing values
#' @param kernelList list of kernels
#' @param names name of kernel
#' @param formula formula
#' @param ... other arguments
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @name kernel.method
#' @rdname kernel.method
#' @aliases comb
#' @aliases check.integer
#' @aliases asOneSidedFormula
#' @aliases splitFormula
#' @aliases computes.Kernel
#' @aliases computes.Kernel.interaction
#' @aliases computes.KernelALL
#' @aliases renames.Kernel
#' @aliases objects.Kernel
#'
#' @export



################################################################################
# Combinatory function
################################################################################

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


################################################################################
# Method checking if a number is an integer (used in polynomial kernel)
################################################################################

#' @rdname kernel.method
#' @export
check.integer <- function(N){
  !grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE))
}


################################################################################
# Method helping to split formula Kernel
################################################################################

#' @rdname kernel.method
#' @export
asOneSidedFormula <- function(object)
{
  if ((mode(object) == "call") && (object[[1L]] == "~")) {
    object <- eval(object, envir = parent.frame())
  }
  if (inherits(object, "formula")) {
    if (length(object) != 2L) {
      stop(gettextf("formula '%s' must be of the form '~expr'",
                    deparse(as.vector(object))), domain = NA)
    }
    return(object)
  }
  do.call("~", list(switch(mode(object), name = , numeric = ,
                           call = object, character = as.name(object), expression = object[[1L]],
                           stop(gettextf("'%s' cannot be of mode '%s'", substitute(object),
                                         mode(object)), domain = NA))))
}


################################################################################
# Method splitting formula Kernel
################################################################################

#' @rdname kernel.method
#' @export
splitFormula <- function(form, sep = "/")
{
  if (inherits(form, "formula") || mode(form) == "call" &&
      form[[1]] == as.name("~"))

    return(splitFormula(form[[length(form)]], sep = sep))
  if (mode(form) == "call" && form[[1]] == as.name(sep))

    return(do.call("c", lapply(as.list(form[-1]), splitFormula,
                               sep = sep)))
  if (mode(form) == "(")
    return(splitFormula(form[[2]], sep = sep))
  if (length(form) < 1)
    return(NULL)
  list(asOneSidedFormula(form))
}


################################################################################
# Method computing single kernel
################################################################################

#' @rdname kernel.method
#' @export
computes.Kernel <- function(x, ind, nameKernel, not.missing = NULL)
{
  nb.kernels <- length(grep(":", nameKernel[[ind]], fixed = TRUE)) + 1
  if (nb.kernels == 1) {
    Z <- as.matrix(x[[ind]]$Z[not.missing, ])
    colnames(Z) <- colnames(x[[ind]]$Z)
    row.names(Z) <- not.missing
    if (x[[ind]]$kernel.function != "gram.matrix") {
      if (isTRUE(x[[ind]]$kernel.scale)) {
        Z2 <- scale(Z)
        Z.scale <- Z2
        kernel.mean <- attr(Z2, "scaled:center")
        kernel.sd <- attr(Z2, "scaled:scale")
        kernel.scale <- TRUE
        attr(Z.scale,"scaled:center") <- NULL
        attr(Z.scale,"scaled:scale") <- NULL
      } else {
        Z2 <- Z
        Z.scale <- NA
        kernel.mean <- NA
        kernel.sd <- NA
        kernel.scale <- FALSE
      }
      K <- kernel.matrix(Z = Z2, whichkernel = x[[ind]]$kernel.function, rho = x[[ind]]$rho, gamma = x[[ind]]$gamma, d = x[[ind]]$d) # REVOIR PARAMETRES DES KERNELS
    } else {
      K <- x[[ind]]$K[not.missing, not.missing]
      Z.scale <- NA
      kernel.mean <- NA
      kernel.sd <- NA
      kernel.scale <- NA
    }
  } else {
    K <- x[[ind]]$K
    Z <- x[[ind]]$Z
    Z.scale <- x[[ind]]$Z.scale
    kernel.mean <- x[[ind]]$kernel.mean
    kernel.sd <- x[[ind]]$kernel.sd
    kernel.scale <- x[[ind]]$kernel.scale
  }
  return(list(K = K, Z = Z, kernel.function = x[[ind]]$kernel.function, Z.scale = Z.scale, kernel.mean = kernel.mean, kernel.sd = kernel.sd, kernel.scale = kernel.scale, kernel.formula = x[[ind]]$kernel.formula, rho = x[[ind]]$rho, gamma = x[[ind]]$gamma, d = x[[ind]]$d, free.parameters = x[[ind]]$free.parameters, type.object = x[[ind]]$type.object))
}


################################################################################
# Method computing interaction between kernels
################################################################################


#' @rdname kernel.method
#' @export
computes.Kernel.interaction <- function(x, ind, nameKernel, not.missing = NULL)
{
  nb.kernels <- length(grep(":", nameKernel[[ind]], fixed = TRUE)) + 1
  if (nb.kernels > 1) {
    num.kernel <- which(nameKernel %in% strsplit(nameKernel[[ind]], ":")[[1]])
    null.kernel = 0
    for (j in 1:length(num.kernel)) {
      null.kernel <- null.kernel + is.null(x[[num.kernel[j]]]$K)
    }
    if (null.kernel == 0) {
      interKernel <- x[[num.kernel[1]]]$K[not.missing, not.missing]
      for (j in 2:length(num.kernel)) {
        interKernel <- interKernel * x[[num.kernel[j]]]$K[not.missing, not.missing]
      }
      K <- interKernel
      kernel.mean <- NULL
      kernel.sd <- NULL
      kernel.function <- NULL
      kernel.scale <- NULL
    } else {
      K <- matrix(NA, nrow = length(not.missing), ncol = length(not.missing))
      rownames(K) <- not.missing
      colnames(K) <- not.missing
      kernel.mean <- x[[ind]]$kernel.mean
      kernel.sd <- x[[ind]]$kernel.sd
      kernel.function <- x[[ind]]$kernel.function
      kernel.scale <- x[[ind]]$kernel.scale
    }
  } else {
    K <- x[[ind]]$K
    kernel.mean <- x[[ind]]$kernel.mean
    kernel.sd <- x[[ind]]$kernel.sd
    kernel.function <- x[[ind]]$kernel.function
    kernel.scale <- x[[ind]]$kernel.scale
  }
  return(list(K = K, Z = x[[ind]]$Z, kernel.function = kernel.function, Z.scale = x[[ind]]$Z.scale, kernel.mean = kernel.mean, kernel.sd = kernel.sd, kernel.scale = kernel.scale, kernel.formula = x[[ind]]$kernel.formula, rho = x[[ind]]$rho, gamma = x[[ind]]$gamma, d = x[[ind]]$d, free.parameters = x[[ind]]$free.parameters, type.object = x[[ind]]$type.object))
}


################################################################################
# Method computing all kernels
################################################################################



#' @rdname kernel.method
#' @export
computes.KernelALL <- function(kernelList, not.missing = NULL){
  kernelList2 <- list()
  kernelList3 <- list()
  for (i in 1:length(kernelList$listKernel)) {
    kernelList2[[i]] <- computes.Kernel(kernelList$listKernel, i, kernelList$nameKernel, not.missing = not.missing)
  }
  for (i in 1:length(kernelList$listKernel)) {
    kernelList3[[i]] <- computes.Kernel.interaction(kernelList2, i, kernelList$nameKernel, not.missing = not.missing)
  }
  kernelList$listKernel <- kernelList3
  return(kernelList)
}


################################################################################
# Method renaming a kernel row/columns
################################################################################



#' @rdname kernel.method
#' @export
renames.Kernel <- function(object, names)
{
  if (!is.null(object$K)) {
    rownames(object$K) <- names
    colnames(object$K) <- names
  }
  if (!is.null(object$Z)) {
    rownames(object$Z) <- names
  }
  return(object)
}


################################################################################
# Method indicating if the kernel is from user own function
################################################################################

#' @rdname kernel.method
#' @export
objects.Kernel <- function(formula)
{
  splitFormula <- strsplit(as.character(formula[2]), "Kernel")[[1]][-1]
  splitFormula <- gsub(" ", "", splitFormula)
  splitFormula <- gsub("\\)\\+", "", splitFormula)
  splitFormula <- gsub("\\(", "", splitFormula)
  splitFormula <- gsub("\\)", "", splitFormula)
  nk <- length(splitFormula)
  out <- list()
  inform.user.matrix <- c()
  for (i in 1:nk) {
    splitForX <- strsplit(splitFormula[i], "x=")
    splitForUserMatrix <- strsplit(splitFormula[i], "kernel.function=\"gram.matrix\"")
    if (length(splitForX[[1]]) == 1) {
      splitForComma <- strsplit(splitFormula[i], ",")[[1]]
      out[[i]] <- splitForComma[1]
    } else {
      splitForComma <- strsplit(splitForX[[1]][2], ",")[[1]]
      out <- c(out, splitForComma[1])
    }
    if (length(splitForUserMatrix[[1]]) == 1 && nchar(splitForUserMatrix[[1]]) + nchar("kernel.function=\"gram.matrix\"") != nchar(splitFormula)) {
      inform.user.matrix[i] <- ""
    } else {
      inform.user.matrix[i] <- "gram.matrix"
    }
  }
  return(list(out, inform.user.matrix))
}



