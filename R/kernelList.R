#' @title List of kernel parts included in the kernel semi parametric model
#'
#'
#' @description internal method for listing all kernel parts included in the model
#'
#' @param formula kernel part formula provided in the \code{kspm} function.
#' @param data data provided in the \code{kspm} function.
#' @param names row names of samples as they are evaluated in \code{kspm} function.
#'
#' @importFrom utils combn
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @export




kernel.list <- function(formula, data, names)
{

  ################################################################################
  # CREATE OBJECTS RETURNED BY THE FUNCTION
  ################################################################################

  # List of kernels computed for the model
  listKernel <- list()
  # Names of these kernels
  nameKernel <- c()
  # List of kernels included in the model
  inclKernel <- c()

  # List of expression used into Kernel() function to characterize the kernel
  listExpres <- list()       # List of all expressions
  mergExpres <- listExpres   # This list will regroup expressions characterizing the same Kernel()

  # Iterative number for expression
  a <- 1
  # Iterative number for kernel
  b <- 1

  # Split the formula
  formula.split <- splitFormula(formula, sep = "+")

  for (i in 1:length(formula.split)) {

    ################################################################################
    # SUBFORMULA IS ONE KERNEL
    ################################################################################

    if (formula.split[[i]][[2]][[1]] == "Kernel") {

      # If this kernel was already evaluated using the same expression
      if (list(formula.split[[i]]) %in% listExpres) {
        # We look for its name and include it in the analysis
        inclKernel <- c(inclKernel, which(nameKernel == paste0("Ker", which(sapply(listExpres, function(x) x == formula.split[[i]])))))
      } else {
        # If this kernel was not already evaluated using the same expression
        # Test if this kernel was already evaluated with another expression

        newKernel <- renames.Kernel(eval(formula.split[[i]][[2]], data), names = names) # Compute the new kernel
        exist <- FALSE
        if (length(listKernel) > 0) {
          for (j in 1:length(listKernel)) {
            # Comparison with existing kernels
            if (ifelse(is.null(dim(newKernel$Z)), 0, dim(newKernel$Z)[2]) != 0 && ifelse(is.null(dim(listKernel[[j]]$Z)), 0, dim(listKernel[[j]]$Z)[2]) != 0 && sum(!(colnames(newKernel$Z) %in% colnames(listKernel[[j]]$Z))) == 0 && sum(!(colnames(newKernel$Z) %in% colnames(listKernel[[j]]$Z))) == 0 && sum(!(colnames(listKernel[[j]]$Z) %in% colnames(newKernel$Z))) == 0 && ifelse(is.null(newKernel$kernel.function), "", newKernel$kernel.function) == ifelse(is.null(listKernel[[j]]$kernel.function), "", listKernel[[j]]$kernel.function) && ifelse(is.null(newKernel$rho), -9, newKernel$rho) == ifelse(is.null(listKernel[[j]]$rho), -9, listKernel[[j]]$rho) && ifelse(is.null(newKernel$gamma), -9, newKernel$gamma) == ifelse(is.null(listKernel[[j]]$gamma), -9, listKernel[[j]]$gamma) && ifelse(is.null(newKernel$d), -9, newKernel$d) == ifelse(is.null(listKernel[[j]]$d), -9, listKernel[[j]]$d) && ifelse(is.null(dim(listKernel[[j]]$Z.scale)), 0, dim(listKernel[[j]]$Z.scale)[2]) == ifelse(is.null(dim(newKernel$Z.scale)), 0, dim(newKernel$Z.scale)[2]) && identical(listKernel[[j]]$Z, newKernel$Z)) {
              warning(paste0("Following kernel: ", as.character(formula.split[[i]][2]), " is considered as similar to: ", as.character(listExpres[[which(nameKernel == paste0("Ker", j))]][2])))
              mergExpres[[which(nameKernel == paste0("Ker", j))]] = c(mergExpres[[which(nameKernel == paste0("Ker", j))]], formula.split[[i]])
              inclKernel <- c(inclKernel, which(nameKernel == paste0("Ker", j)))
              exist <- TRUE
            }
          }
        }
        # If this kernel was not already evaluated using a different expression
        # We added the new kernel in the lists
        if (exist == FALSE) {
          # List of expressions
          listExpres[[a]] <- formula.split[[i]]
          mergExpres[[a]] <- formula.split[[i]]
          # Lists of kernels
          listKernel[[b]] <- newKernel
          inclKernel <- c(inclKernel, b)
          nameKernel <- c(nameKernel, paste0("Ker", a))
          # Modify iterative numbers
          a <- a + 1
          b <- b + 1
        }
      }
    }

    ################################################################################
    # SUBFORMULA IS A COMPLEXE INTERACTION OF KERNELS
    ################################################################################

    if (formula.split[[i]][[2]][[1]] == "*") {
      # The formula is splitted into subformulae
      subSplitFormula <- splitFormula(formula.split[[i]], sep = "*")
      # Search for all kernels included in formula.split[[i]]
      num <- c()   # Kernels included in formula.split[[i]]
      for (j in 1:length(subSplitFormula)) {
        # If this kernel was already evaluated using the same expression
        if (list(subSplitFormula[[j]]) %in% listExpres) {
          # We look for its name and include it in the analysis
          num <- c(num, which(sapply(listExpres, function(x) x == subSplitFormula[[j]])))
          inclKernel <- c(inclKernel, which(sapply(listExpres, function(x) x == subSplitFormula[[j]])))
        } else {
          # If this kernel was not already evaluated using the same expression
          # Test if this kernel was already evaluated with another expression
          newKernel <- renames.Kernel(eval(subSplitFormula[[j]][[2]], data), names = names)
          exist <- FALSE
          if (length(listKernel) > 0) {
            for (l in 1:length(listKernel)) {
              # Comparison with existing kernels
              if (ifelse(is.null(dim(newKernel$Z)), 0, dim(newKernel$Z)[2]) != 0 && ifelse(is.null(dim(listKernel[[l]]$Z)), 0, dim(listKernel[[l]]$Z)[2]) != 0 && sum(!(colnames(newKernel$Z) %in% colnames(listKernel[[l]]$Z))) == 0 && sum(!(colnames(newKernel$Z) %in% colnames(listKernel[[l]]$Z))) == 0 && sum(!(colnames(listKernel[[l]]$Z) %in% colnames(newKernel$Z))) == 0 && ifelse(is.null(newKernel$kernel.function), "", newKernel$kernel.function) == ifelse(is.null(listKernel[[l]]$kernel.function), "", listKernel[[l]]$kernel.function) && ifelse(is.null(newKernel$rho), -9, newKernel$rho) == ifelse(is.null(listKernel[[l]]$rho), -9, listKernel[[l]]$rho) && ifelse(is.null(newKernel$gamma), -9, newKernel$gamma) == ifelse(is.null(listKernel[[l]]$gamma), -9, listKernel[[l]]$gamma) && ifelse(is.null(newKernel$d), -9, newKernel$d) == ifelse(is.null(listKernel[[l]]$d), -9, listKernel[[l]]$d) && ifelse(is.null(dim(listKernel[[l]]$Z.scale)), 0, dim(listKernel[[l]]$Z.scale)[2]) == ifelse(is.null(dim(newKernel$Z.scale)), 0, dim(newKernel$Z.scale)[2]) && identical(listKernel[[l]]$K, newKernel$K)) {
                num <- c(num, l)
                mergExpres[[which(nameKernel == paste0("Ker", l))]] = c(mergExpres[[which(nameKernel == paste0("Ker", l))]], subSplitFormula[[j]])
                inclKernel <- c(inclKernel, which(nameKernel == paste0("Ker", l)))
                exist <- TRUE
              }
            }
          }
          # If this kernel was not already evaluated using a different expression
          # We added the new kernel in the lists
          if (exist == FALSE) {
            # List of expressions
            listExpres[[a]] <- subSplitFormula[[j]]
            mergExpres[[a]] <- subSplitFormula[[j]]
            # List of kernels
            listKernel[[b]] <- renames.Kernel(eval(subSplitFormula[[j]][[2]], data), names = names)
            inclKernel <- c(inclKernel, b)
            nameKernel <- c(nameKernel, paste0("Ker", a))
            num <- c(num, a)
            # Modify iterative numbers
            b <- b + 1
            a <- a + 1
          }
        }
      }
      # Compute all possible interactions with kernels included in formula.split[[i]]
      numU <- unique(num)
      for (j in 2:length(numU)) {   # j represents the number of kernels in the interaction
        X <- combn(length(numU), j) # X represents the combination of j kernels
        for (l in 1:(dim(X)[2])) {  # l represents one of these combination
          # Name of the interaction
          namelistKernel <- paste0("Ker", sort(numU[X[,l]]), collapse = ":")
          if (!(namelistKernel %in% nameKernel)) {
            nameKernel <- c(nameKernel, namelistKernel)
            listKernel[[b]] <- list(K = NULL, Z = NULL, kernel.function = NULL, Z.scale = NULL, kernel.mean = NULL, kernel.sd = NULL, kernel.scale = NULL, kernel.formula = NULL, rho = NULL, gamma = NULL, d = NULL)
            inclKernel <- c( inclKernel, b)
            b <- b + 1
          }
        }
      }
    }

    ################################################################################
    # SUBFORMULA IS A SIMPLE INTERACTION OF KERNELS
    ################################################################################

    if (formula.split[[i]][[2]][[1]] == ":") {
      # The formula is splitted into subformulae
      subSplitFormula <- splitFormula(formula.split[[i]], sep = ":")
      # Search for all kernels included in formula.split[[i]]
      num <- c()   # Kernels included in formula.split[[i]]
      for (j in 1:length(subSplitFormula)) {
        # If this kernel was already evaluated using the same expression
        if (list(subSplitFormula[[j]]) %in% listExpres) {
          num <- c(num, which(sapply(listExpres, function(x) x == subSplitFormula[[j]])))
        } else {
          # If this kernel was not already evaluated using the same expression
          # Test if this kernel was already evaluated with another expression
          # newKernel <- renames.Kernel(eval(subSplitFormula[[j]][[2]], data, envir = parent.frame(2)), names = names)
          newKernel <- renames.Kernel(eval(subSplitFormula[[j]][[2]], data), names = names)
          exist <- FALSE
          if (length(listKernel) > 0) {
            for (l in 1:length(listKernel)) {
              # Comparison with existing kernels
              if (ifelse(is.null(dim(newKernel$Z)), 0, dim(newKernel$Z)[2]) != 0 && ifelse(is.null(dim(listKernel[[l]]$Z)), 0, dim(listKernel[[l]]$Z)[2]) != 0 && sum(!(colnames(newKernel$Z) %in% colnames(listKernel[[l]]$Z))) == 0 && sum(!(colnames(newKernel$Z) %in% colnames(listKernel[[l]]$Z))) == 0 && sum(!(colnames(listKernel[[l]]$Z) %in% colnames(newKernel$Z))) == 0 && ifelse(is.null(newKernel$kernel.function), "", newKernel$kernel.function) == ifelse(is.null(listKernel[[l]]$kernel.function), "", listKernel[[l]]$kernel.function) && ifelse(is.null(newKernel$rho), -9, newKernel$rho) == ifelse(is.null(listKernel[[l]]$rho), -9, listKernel[[l]]$rho) && ifelse(is.null(newKernel$gamma), -9, newKernel$gamma) == ifelse(is.null(listKernel[[l]]$gamma), -9, listKernel[[l]]$gamma) && ifelse(is.null(newKernel$d), -9, newKernel$d) == ifelse(is.null(listKernel[[l]]$d), -9, listKernel[[l]]$d) && ifelse(is.null(dim(listKernel[[l]]$Z.scale)), 0, dim(listKernel[[l]]$Z.scale)[2]) == ifelse(is.null(dim(newKernel$Z.scale)), 0, dim(newKernel$Z.scale)[2]) && identical(listKernel[[l]]$K, newKernel$K)) {
                inclKernel <- c(inclKernel, which(nameKernel == paste0("Ker", l)))
                warning(paste0("Following kernel: ", as.character(subSplitFormula[[j]][2]), " is considered as similar to: ", as.character(listExpres[[which(nameKernel == paste0("Ker", l))]][2])))
                mergExpres[[which(nameKernel == paste0("Ker", l))]] = c(mergExpres[[which(nameKernel == paste0("Ker", l))]], subSplitFormula[[j]])
                exist <- TRUE
              }
            }
          }
          # If this kernel was not already evaluated using a different expression
          # We added the new kernel in the lists
          if (exist == FALSE) {
            # List of expressions
            listExpres[[a]] <- subSplitFormula[[j]]
            mergExpres[[a]] <- subSplitFormula[[j]]

            # List of kernels
            # listKernel[[b]] <- renames.Kernel(eval(subSplitFormula[[j]][[2]], data, envir = parent.frame(2)), names = names)
            listKernel[[b]] <- renames.Kernel(eval(subSplitFormula[[j]][[2]], data), names = names)

            nameKernel <- c(nameKernel, paste0("Ker", a))
            num <- c(num, a)
            # Modify iterative numbers
            b <- b + 1
            a <- a + 1
          }
        }
      }
      # Compute the interaction with kernels included in formula.split[[i]]
      numU <- unique(num)
      X <- combn(length(numU), length(numU)) # X represents the combination of the kernels
      # Name of the interaction
      namelistKernel <- paste0("Ker", sort(numU[X[,1]]), collapse = ":")
      # If this interaction was not already computed
      # We create the corresponding interaction kernel
      if (!(namelistKernel %in% nameKernel)) {
        nameKernel <- c(nameKernel, namelistKernel)
        listKernel[[b]] <- list(K = NULL, Z = NULL, kernel.function = NULL, Z.scale = NULL, kernel.mean = NULL, kernel.sd = NULL, kernel.scale = NULL, kernel.formula = NULL, rho = NULL, gamma = NULL, d = NULL)
        inclKernel <- c( inclKernel, b)
        b <- b + 1
      }
    }


    ################################################################################
    # SUBFORMULA IS A LIMITED INTERACTION OF KERNELS
    ################################################################################

    if (formula.split[[i]][[2]][[1]] == "^") {
      # The formula is splitted into subformulae
      insideFormula <- splitFormula(formula.split[[i]], sep = "^")
      # Stop if "*" or ":" are founded inside exposured term
      if (length(grep("*", as.character(insideFormula[[1]][2]), fixed = TRUE)) + length(grep(":", as.character(insideFormula[[1]][2]), fixed = TRUE)) > 0) {
        stop("'*' and ':' are not allowed into exposured term in kernel formula")
      }
      subSplitFormula <- splitFormula(insideFormula[[1]], sep = "+")
      exposureTerm <- insideFormula[[2]][[2]]
      # Search for all kernels included in formula.split[[i]]
      num <- c()   # Kernels included in formula.split[[i]]
      for (j in 1:length(subSplitFormula)) {
        # If this kernel was already evaluated using the same expression
        if (list(subSplitFormula[[j]]) %in% listExpres) {
          # We look for its name and include it in the analysis
          num <- c(num, which(sapply(listExpres, function(x) x == subSplitFormula[[j]])))
          inclKernel <- c(inclKernel, which(sapply(listExpres, function(x) x == subSplitFormula[[j]])))
        } else {
          # If this kernel was not already evaluated using the same expression
          # Test if this kernel was already evaluated with another expression
          # newKernel <- renames.Kernel(eval(subSplitFormula[[j]][[2]], data, envir = parent.frame(2)), names = names)
          newKernel <- renames.Kernel(eval(subSplitFormula[[j]][[2]], data), names = names)

          exist <- FALSE
          if (length(listKernel) > 0) {
            for (l in 1:length(listKernel)) {
              # Comparison with existing kernels
              if (ifelse(is.null(dim(newKernel$Z)), 0, dim(newKernel$Z)[2]) != 0 && ifelse(is.null(dim(listKernel[[l]]$Z)), 0, dim(listKernel[[l]]$Z)[2]) != 0 && sum(!(colnames(newKernel$Z) %in% colnames(listKernel[[l]]$Z))) == 0 && sum(!(colnames(newKernel$Z) %in% colnames(listKernel[[l]]$Z))) == 0 && sum(!(colnames(listKernel[[l]]$Z) %in% colnames(newKernel$Z))) == 0 && ifelse(is.null(newKernel$kernel.function), "", newKernel$kernel.function) == ifelse(is.null(listKernel[[l]]$kernel.function), "", listKernel[[l]]$kernel.function) && ifelse(is.null(newKernel$rho), -9, newKernel$rho) == ifelse(is.null(listKernel[[l]]$rho), -9, listKernel[[l]]$rho) && ifelse(is.null(newKernel$gamma), -9, newKernel$gamma) == ifelse(is.null(listKernel[[l]]$gamma), -9, listKernel[[l]]$gamma) && ifelse(is.null(newKernel$d), -9, newKernel$d) == ifelse(is.null(listKernel[[l]]$d), -9, listKernel[[l]]$d) && ifelse(is.null(dim(listKernel[[l]]$Z.scale)), 0, dim(listKernel[[l]]$Z.scale)[2]) == ifelse(is.null(dim(newKernel$Z.scale)), 0, dim(newKernel$Z.scale)[2]) && identical(listKernel[[l]]$K, newKernel$K)) {
                num <- c(num, l)
                warning(paste0("Following kernel: ", as.character(subSplitFormula[[j]][2]), " is considered as similar to: ", as.character(listExpres[[which(nameKernel == paste0("Ker", l))]][2])))
                mergExpres[[which(nameKernel == paste0("Ker", l))]] = c(mergExpres[[which(nameKernel == paste0("Ker", l))]], subSplitFormula[[j]])
                inclKernel <- c(inclKernel, which(nameKernel == paste0("Ker", l)))
                exist <- TRUE
              }
            }
          }
          # If this kernel was not already evaluated using a different expression
          # We added the new kernel in the lists
          if (exist == FALSE) {
            # List of expressions
            listExpres[[a]] <- subSplitFormula[[j]]
            mergExpres[[a]] <- subSplitFormula[[j]]
            # List of kernels
            # listKernel[[b]] <- renames.Kernel(eval(subSplitFormula[[j]][[2]], data, envir = parent.frame(2)), names = names)
            listKernel[[b]] <- renames.Kernel(eval(subSplitFormula[[j]][[2]], data), names = names)

            inclKernel <- c( inclKernel, b)
            nameKernel <- c(nameKernel, paste0("Ker", a))
            num <- c(num, a)
            # Modify iterative numbers
            b <- b + 1
            a <- a + 1
          }
        }
      }
      # Compute all possible interactions with kernels included in formula.split[[i]]
      numU <- unique(num)
      # Test if exposureTerm is plausible: if not => Stop
      if (exposureTerm > length(numU)) {
        stop("exposure term in kernel formula could not exceed the number of kernels in the scope")
      }
      for (j in 2:exposureTerm) {   # j represents the number of kernels in the interaction
        X <- combn(length(numU), j) # X represents the combination of j kernels
        for (l in 1:(dim(X)[2])) {  # l represents one of these combination
          # Name of the interaction
          namelistKernel <- paste0("Ker", sort(numU[X[,l]]), collapse = ":")
          # If this interaction was not already computed
          # We create the corresponding interaction kernel
          if (!(namelistKernel %in% nameKernel)) {
            nameKernel <- c(nameKernel, namelistKernel)
            listKernel[[b]] <- list(K = NULL, Z = NULL, kernel.function = NULL, Z.scale = NULL, kernel.mean = NULL, kernel.sd = NULL, kernel.scale = NULL, kernel.formula = NULL, rho = NULL, gamma = NULL, d = NULL)
            inclKernel <- c( inclKernel, b)
            b <- b + 1
          }
        }
      }
    }
  }


  ################################################################################
  # COUTPUT
  ################################################################################

  # unique kernels
  singKernel1 <- which(unique(nameKernel) %in% unique(nameKernel)[-(grep(":", nameKernel))])
  if (length(singKernel1) == 0 ) {
    singKernel <- unique(inclKernel)
  } else {
    singKernel <- singKernel1
  }

  # constraints for free parameters
  free.parameters.contraint <- c()
  for (i in 1:length(listKernel)) {
    free.parameters <- c()
    if (!is.null(listKernel[[i]]$kernel.function)) {
      if (listKernel[[i]]$kernel.function == "gaussian" & is.null(listKernel[[i]]$rho)) {
        free.parameters <- c(free.parameters, "rho")
        free.parameters.contraint <- c(free.parameters.contraint, 0)
      }
      if (listKernel[[i]]$kernel.function == "polynomial" & is.null(listKernel[[i]]$rho)) {
        free.parameters <- c(free.parameters, "rho")
        free.parameters.contraint <- c(free.parameters.contraint, 0)
      }
      if (listKernel[[i]]$kernel.function == "sigmoid" & is.null(listKernel[[i]]$rho)) {
        free.parameters <- c(free.parameters, "rho")
        free.parameters.contraint <- c(free.parameters.contraint, -Inf)
      }
      if (listKernel[[i]]$kernel.function == "polynomial" & is.null(listKernel[[i]]$gamma)) {
        free.parameters <- c(free.parameters, "gamma")
        free.parameters.contraint <- c(free.parameters.contraint, -Inf)
      }
      if (listKernel[[i]]$kernel.function == "sigmoid" & is.null(listKernel[[i]]$gamma)) {
        free.parameters <- c(free.parameters, "gamma")
        free.parameters.contraint <- c(free.parameters.contraint, -Inf)
      }
      if (listKernel[[i]]$kernel.function == "inverse.quadratic" & is.null(listKernel[[i]]$gamma)) {
        free.parameters <- c(free.parameters, "gamma")
        free.parameters.contraint <- c(free.parameters.contraint, 0)
      }
    }
    if (length(free.parameters) == 0) {
      listKernel[[i]]$free.parameters <- ""
    } else {
      listKernel[[i]]$free.parameters <- free.parameters
    }
  }

  return(list(listKernel = listKernel, nameKernel = nameKernel, singKernel = singKernel, inclKernel = unique(inclKernel), listExpres = listExpres, mergExpres = mergExpres, free.parameters.contraint = free.parameters.contraint))
}

