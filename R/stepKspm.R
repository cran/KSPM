#' @title Choose a model by AIC or BIC in a Stepwise Algorithm
#'
#' @description Performs stepwise model selection for Kernel Semi Parametric Model by AIC or BIC.
#' @param object an object of class "kspm" with only one kernel.
#' @param data data.
#' @param linear.lower one side formula corresponding to the smallest set of variables that should be included in the linear part of the model.
#' @param linear.upper one side formula corresponding to the largest set of variables that may be included in the linear part of the model.
#' @param kernel.lower one side formula corresponding to the smallest set of variables that should be included in the kernel part of the model.
#' @param kernel.upper one side formula corresponding to the  largest set of variables that may be included in the kernel part of the model.
#' @param direction the mode of stepwise search, can be one of "both" (default), "backward", or "forward".
#' @param k type of information criteria used for the variable selection. If \code{k=2} AIC is used (default), if \code{k=log(n)}, BIC is used instead.
#' @param kernel.param define if hyperparameters should be fixed (\code{"fixed"}) or reestimated at each iteration (\code{"change"}). Tu use the last option, hyperparameter of model provided in \code{object} should have been estimated by the model.
#' @param trace integer. If positive, information is printed during the running of step.kspm. Larger values may give more information on the fitting process.
#'
#'
#' @return \code{stepKSPM} returns the selected model.
#'
#'
#' @details This procedure may be done on \code{kspm} object defined with only one kernel part and for which a data frame including all variables was provided. Selection may be done on linear part only, on kernel part only or on both at the same time. To perform selection on linear (resp. kernel) part only, \code{kernel.lower} and \code{kernel.upper} (resp. \code{linear.lower} and \code{linear.upper}) should contain all the variables that should stay in the model for kernel (resp. linear) part.
#'
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @importFrom stats update.formula
#'
#'
#' @seealso \link{extractAIC.kspm}
#'
#' @examples
#' x <- 1:15
#' z1 <- runif(15, 1, 6)
#' z2 <- rnorm(15, 1, 4)
#' z3 <- rnorm(15, 6, 2)
#' z4 <- runif(15, -10, 2)
#' y <- 3*x + (z1 + z2)^2 + rnorm(15, 0, 2)
#' dfrm <- data.frame(x = x, z1 = z1, z2 = z2, z3 = z3, z4 = z4, y = y)
#' fit <- kspm(y, linear = ~ x, kernel = ~ Kernel(~ z1 + z2 + z3 + z4,
#' kernel.function = "polynomial", d= 2, rho = 1, gamma = 0), data = dfrm)
#' stepKSPM(fit, k = 2, data = dfrm)
#'
#' @export




stepKSPM <- function(object, data = NULL, linear.lower = NULL, linear.upper = NULL, kernel.lower = NULL, kernel.upper = NULL, direction = "both", k = 2, kernel.param = "fixed", trace = TRUE){

  ################################################################################
  # CALL
  ################################################################################

  call.step <- match.call(expand.dots = FALSE)
  call.names.step <- paste0(call.step)
  call.order.step <- match(c("linear.lower", "linear.upper", "kernel.lower", "kernel.upper", "data"), names(call.step), 0L)
  call.order <- match(c("linear", "kernel", "response"), names(object$call), 0L)

  ################################################################################
  # SOME CHECKS FOR ARGUMENTS
  ################################################################################

  if (!inherits(object, "kspm")) {
    stop("object should be of class 'kspm'.")
  }
  if (!is.null(linear.lower) & !inherits(linear.lower, "formula")) {
    stop("linear.lower should be a formula.")
  }
  if (!is.null(linear.upper) & !inherits(linear.upper, "formula")) {
    stop("linear.upper should be a formula.")
  }
  if (!is.null(kernel.lower) & !inherits(kernel.lower, "formula")) {
    stop("kernel.lower should be a formula.")
  }
  if (!is.null(kernel.upper) & !inherits(kernel.upper, "formula")) {
    stop("kernel.upper should be a formula.")
  }
  if (!(is.data.frame(data) | is.matrix(data))) {
    stop(paste0(data, " is not a data.frame or a matrix."))
  }
  if (!(trace %in% c(TRUE, FALSE))) {
    stop(paste0(trace, " is not a valid argument for 'trace' option. Choose TRUE or FALSE."))
  }
  if (!(direction %in% c("forward", "backward", "both"))) {
    stop(paste0(direction, " is not a valid argument for 'direction' option. Choose 'forward', 'backward' or 'both'."))
  }
  if (!(kernel.param %in% c("fixed", "change"))) {
    stop(paste0(kernel.param, " is not a valid argument for 'kernel.param' option, Choose 'fixed' or 'change'."))
  }

  ################################################################################
  # CASES FOR WHICH stepKSPM HAS NOT BEEN IMPLEMENTED YET
  ################################################################################

  # if the kernel object has been defined using Gram matrix (user kernel function), the stepKSPM function could not be used
  informUserMatrix <- objects.Kernel(object$call[call.order[2]][[1]])[[2]]
  if (informUserMatrix == "user.matrix") {
    stop("step could not worked when user kernel matrix is provided")
  }
  # if there are more than 1 kernel part in the object, the stepKSPM function could not be used
  if (length(object$lambda) > 1) {
    stop("step could not worked for more than 1 kernel part")
  }

  ################################################################################
  # CONTROL PARAMETERS
  ################################################################################

  # control parameters are those used for computing initial model
  kernel.control <- object$control

  ################################################################################
  # INITIALISATION LINEAR PART
  ################################################################################

  # lower formula
  # if NULL, lower is intercept only, else is linear.lower
  if (is.null(linear.lower)) {
    linear.lower <- ~1
  }

  # upper formula
  # if NULL, upper is linear formula from the object, else is linear.upper
  if (is.null(linear.upper)) {
    linear.upper <- object$linear.formula
  } else {# if is non NULL but X has been defined by a matrix in object, we use columns of X instead
    if (is.matrix(eval(object$call[call.order[1]][[1]]))) {
      warning(paste0("linear part in initial model has been defined by a matrix: terms in linear.upper are ignored and replaced by the elements of the matrix ", object$call[call.order[1]][[1]]))
      linear.upper <- object$linear.formula
    }
  }

  # design matrix for candidate variables for linear part
  # if linear part of object has been defined by a matrix, we used it and removed the intercept if needed
  if (is.matrix(eval(object$call[call.order[1]][[1]]))) {
    if (colnames(object$X)[1] == "(Intercept)") {
      data.linear <- data.frame(object$X[, -1])
    } else {
      data.linear <- data.frame(object$X)
    }
  } else {# else we use the formula
    if (is.null(data)) {
      data.linear <- model.frame(linear.upper)
    } else {
      data.linear <- model.frame(linear.upper, data)
    }
  }

  # names of linear variables
  linear.current.names <- colnames(object$X)[colnames(object$X) != "(Intercept)"]
  linear.lower.names <- names(model.frame(linear.lower, data.linear))
  linear.upper.names <- names(data.linear)

  # some checks: coherence between current model and lower/upper fomulae
  # linear.lower.names in linear.current.names
  if (sum(!(linear.lower.names %in% linear.current.names)) > 0) {
    stop("the linear.lower formula has a term non included in the model")
  }
  # linear.current.names in linear.current.names
  if (sum(!(linear.lower.names %in% linear.upper.names)) > 0) {
    stop("the model has a term non included in the linear.upper formula")
  }

  ################################################################################
  # INITIALISATION KERNEL PART
  ################################################################################

  # useful information
  kernelObject <- objects.Kernel(object$call[call.order[2]][[1]])[[1]]

  # lower formula
  # if NULL, lower is intercept only, else is kernel.lower
  if (is.null(kernel.lower)) {
    kernel.lower <- ~1
  }

  # upper formula
  # if NULL, upper is kernel formula from the object, else is kernel.upper
  if (is.null(kernel.upper)) {
    kernel.upper <- object$kernel.info$Ker1$kernel.formula
  } else {# if is non NULL but Z has been defined by a matrix in object, we use columns of Z instead
    if (is.matrix(eval(parse(text = kernelObject[[1]]), envir = parent.frame()))) {
      warning(paste0("kernel part in initial model has been defined by a matrix: terms in kernel.upper are ignored and replaced by the elements of the matrix ", kernelObject[[1]]))
      linear.upper <- object$kernel.info$Ker1$kernel.formula
    }
  }

  # design matrix for candidate variables for kernel part
  # if kernel part of object has been defined by a matrix, we used it
  if (is.matrix(eval(parse(text = kernelObject[[1]]), envir = parent.frame())))  {
    data.kernel <- data.frame(object$kernel.info$Ker1$Z)
  } else {# else we use the formula
    if (is.null(data)) {
      data.kernel <- model.frame(kernel.upper)
    } else {
      data.kernel <- model.frame(kernel.upper, data)
    }
  }

  # names of kernel variables
  kernel.current.names <- colnames(object$kernel.info$Ker1$Z)
  kernel.lower.names <- names(model.frame(kernel.lower, data.kernel))
  kernel.upper.names <- names(data.kernel)

  # some checks: coherence between current model and lower/upper fomulae
  # kernel.lower.names in kernel.current.names
  if (sum(!(kernel.lower.names %in% kernel.current.names)) > 0) {
    stop("the kernel.lower formula has a term non included in the model")
  }
  # kernel.current.names in kernel.current.names
  if (sum(!(kernel.current.names %in% kernel.upper.names)) > 0) {
    stop("the model has a term non included in the kernel.upper formula")
  }

  #-------------------------------------------------------------------------------
  # DEFINE THE KERNEL OPTIONS/PARAMETERS
  #-------------------------------------------------------------------------------

  # kernel options that cannot be changed during the step procedure
  kernel.function <- object$kernel.info$Ker1$kernel.function
  kernel.scale <- object$kernel.info$Ker1$kernel.scale
  kernel.d <- object$kernel.info$Ker1$d

  # Other fixed parameters if option " kernel.param="fixed" "
  if (kernel.param == "fixed") {
    kernel.rho <- object$kernel.info$Ker1$rho
    kernel.gamma <- object$kernel.info$Ker1$gamma
  } else if (kernel.param == "change") { # Parameters that may be modified during the step procedure
    # rho option
    if ("rho" %in% object$kernel.info$Ker1$free.parameter) {
      kernel.rho <- NULL
    } else {
      kernel.rho <- object$kernel.list$Ker1$rho[1]
      if (!is.null(kernel.rho)) {
        warning(paste0("kernel.param = 'change' is ignored because initial model was computed with fixed rho. rho = ", kernel.rho, " instead."))
      }
    }
    # gamma option
    if ("gamma" %in% object$kernel.info$Ker1$free.parameter) {
      kernel.gamma <- NULL
    } else {
      kernel.gamma <- object$kernel.list$Ker1$gamma[1]
      if (!is.null(kernel.gamma)) {
        warning(paste0("kernel.param = 'change' is ignored because initial model was computed with fixed gamma. gamma = ", kernel.gamma, "instead."))
      }
    }
  }

  ################################################################################
  # CHECK MISSING DATA
  ################################################################################

  # data
  dfrm <- na.omit(cbind(data.linear, data.kernel))

  # data should have the same lenth during all the step procedure
  if (dim(dfrm)[1] < object$n) {
    stop("the number of rows has been changed: remove missing data ?")
  }

  ################################################################################
  # STEP PROCEDURE
  ################################################################################

  #-------------------------------------------------------------------------------
  # INITIALISATION
  #-------------------------------------------------------------------------------

  # initial model
  m0 <- object

  # initial AIC
  aic <- extractAIC.kspm(m0, k = k)

  # print initial information (if asked for)
  if (trace) {
    cat(paste0("Start: AIC = ", round(aic,1)))
    cat(paste0("\nLinear part: ~ ", m0$linear.formula)[2])
    cat(paste0("\nKernel part: ~ ", m0$kernel.info$Ker1$kernel.formula)[2])
    cat("\n\n")
  }

  # initialisation of loop parameters
  step <- 0
  continue <- TRUE

  # loop
  while (continue == TRUE) {

    #-------------------------------------------------------------------------------
    # DEFINITION OF VARIABLES TO ADD OR TO REMOVE
    #-------------------------------------------------------------------------------

    # current variable(s)
    linear.current <- names(model.frame(m0$linear.formula, dfrm))
    kernel.current <- NULL
    if (!is.null(m0$kernel.info)) {
      kernel.current <- names(model.frame(m0$kernel.info$Ker1$kernel.formula, dfrm))
    }

    # variable(s) to add
    if (direction == "backward") {
      # if backward, no variable to add
      linear.add <- c()
      kernel.add <- c()
    } else {
      # if forward or both, potential variable(s) to add
      linear.add <- linear.upper.names[!(linear.upper.names %in% linear.current)]
      kernel.add <- kernel.upper.names[!(kernel.upper.names %in% kernel.current)]
    }

    # variable(s) to remove
    if (direction == "forward") {
      # if forward, no variable to remove
      linear.remove <- c()
      kernel.remove <- c()
    } else {
      # if backward or both, potential variable(s) to remove
      linear.remove <- linear.current[!(linear.current %in% linear.lower.names)]
      kernel.remove <- kernel.current[!(kernel.current %in% kernel.lower.names)]
    }

    # if nothing to add or remove the procedure stops, else continues
    if (length(linear.add) + length(linear.remove) + length(kernel.add) + length(kernel.remove) == 0) {
      continue <- FALSE
    } else {

      #-------------------------------------------------------------------------------
      # COMPUTATION OF POSSIBLE MODELS
      #-------------------------------------------------------------------------------

      ### list of models that adds a variable in the linear part
      model.add.linear <- list()

      # if no variable to add into the linear part
      if (length(linear.add) == 0) {
        aic.linear.add <- data.frame(Part = character(), AIC = double())
      } else {# if at least one variable to add into the linear part
        # built the table
        aic.linear.add <- data.frame(Part = rep("linear", length(linear.add)), AIC = rep(NA, length(linear.add)))
        # rownames of the table are candidate variables
        row.names(aic.linear.add) <- paste0("+ ", linear.add)
        # compute model and AIC for each candidate variable
        for (i in 1:length(linear.add)) {
          # new linear formula
          formula.new <- update.formula(old = m0$linear.formula, new = as.formula(paste0("~ . + 1 +", linear.add[i])))
          # new model
          fn <- as.formula(eval(paste0("~Kernel(",paste(m0$kernel.info$Ker1$kernel.formula, collapse = ""),", scale = ", kernel.scale, ", kernel.function = '", kernel.function, "', rho = ", kernel.rho, ", gamma = ", kernel.gamma, ", d = ", kernel.d, ")")))
          m.new <- kspm(response = object$Y, linear = formula.new, kernel = fn , data = dfrm, level = 0, control = kernel.control)
          model.add.linear[[i]] <- m.new
          # AIC
          aic.linear.add[i, "AIC"] <- extractAIC.kspm(m.new, k = k)
        }
      }


      ### list of models that adds a variable in the kernel part
      model.add.kernel <- list()

      # if no variable to add into the kernel part
      if (length(kernel.add) == 0) {
        aic.kernel.add <- data.frame(Part = character(), AIC = double())
      } else {# if at least one variable to add into the kernel part
        # built the table
        aic.kernel.add <- data.frame(Part = rep("kernel", length(kernel.add)), AIC = rep(NA, length(kernel.add)))
        # rownames of the table are candidate variables
        row.names(aic.kernel.add) <- paste0("+ ", kernel.add)
        # compute model and AIC for each candidate variable
        for (i in 1:length(kernel.add)) {
          # new kernel formula
          if (!is.null(m0$kernel.info)) { # if current model has a kernel part
            formula.new <- update.formula(old = m0$kernel.info$Ker1$kernel.formula, new = as.formula(paste0("~ . + 1 +", kernel.add[i])))
          } else {# if current model has not a kernel part
            formula.new <- as.formula(paste0("~  1 +", kernel.add[i]))
          }
          # new model
          fn <- as.formula(eval(paste0("~Kernel(",paste(formula.new, collapse = ""), ", scale = ", kernel.scale, ", kernel.function = '", kernel.function, "', rho = ", kernel.rho, ", gamma = ", kernel.gamma, ", d = ", kernel.d, ")")))
          m.new <- kspm(response = object$Y, linear = m0$linear.formula, kernel = fn , data = dfrm, level = 0, control = kernel.control)
          model.add.kernel[[i]] <- m.new
          # AIC
          aic.kernel.add[i, "AIC"] <- extractAIC.kspm(m.new, k = k)
        }
      }


      ### list of models that removes a variable in the linear part
      model.remove.linear <- list()

      # if no variable to remove from the linear part
      if (length(linear.remove) == 0) {
        aic.linear.remove <- data.frame(Part = character(), AIC = double())
      } else {# if at least one variable to remove from the linear part
        # built the table
        aic.linear.remove <- data.frame(Part = rep("linear", length(linear.remove)), AIC = rep(NA, length(linear.remove)))
        # rownames of the table are candidate variables
        row.names(aic.linear.remove) <- paste0("- ", linear.remove)
        # compute model and AIC for each candidate variable
        for (i in 1:length(linear.remove)) {
          # new linear formula
          formula.new <- update.formula(old = m0$linear.formula, new = as.formula(paste0("~ . - ", linear.remove[i])))
          # new model
          fn <- as.formula(eval(paste0("~Kernel(",paste(m0$kernel.info$Ker1$kernel.formula, collapse = ""),", scale = ", kernel.scale, ", kernel.function = '", kernel.function, "', rho = ", kernel.rho, ", gamma = ", kernel.gamma, ", d = ", kernel.d, ")")))
          m.new <- kspm(response = object$Y, linear = formula.new, kernel = fn , data = dfrm, level = 0, control = kernel.control)
          model.remove.linear[[i]] <- m.new
          # AIC
          aic.linear.remove[i, "AIC"] <- extractAIC.kspm(m.new, k = k)
        }
      }


      ### list of models that removes a variable in the kernel part
      model.remove.kernel <- list()

      # if no variable to remove from the kernel part
      if (length(kernel.remove) == 0) {
        aic.kernel.remove <- data.frame(Part = character(), AIC = double())
      } else if (length(kernel.remove) == 1) { # if only one candidate variable
        # built the table
        aic.kernel.remove <- data.frame(Part = rep("kernel", length(kernel.remove)), AIC = rep(NA, length(kernel.remove)))
        # rownames of the table are candidate variables
        row.names(aic.kernel.remove) <- paste0("- ", kernel.remove)
        # if new model has a NULL kernel
        if (length(kernel.lower.names) == 0) {
          # new model
          m.new <- kspm(response = object$Y, linear = m0$linear.formula, kernel = NULL , data = dfrm, level = 0, control = kernel.control)
          model.remove.kernel[[1]] <- m.new
          # AIC
          aic.kernel.remove[1, "AIC"] <- extractAIC.kspm(m.new, k = k)
        } else {# if new model has not a NULL kernel
          # new kernel formula
          formula.new <- update.formula(old = m0$kernel.info$Ker1$kernel.formula, new = as.formula(paste0("~ . - ", kernel.remove[1])))
          # new model
          fn <- as.formula(eval(paste0("~Kernel(",paste(formula.new, collapse = ""),", scale = ", kernel.scale, ", kernel.function = '", kernel.function, "', rho = ", kernel.rho, ", gamma = ", kernel.gamma, ", d = ", kernel.d, ")")))
          m.new <- kspm(response = object$Y, linear = m0$linear.formula, kernel = fn , data = dfrm, level = 0, control = kernel.control)
          model.remove.kernel[[1]] <- m.new
          # AIC
          aic.kernel.remove[1, "AIC"] <- extractAIC.kspm(m.new, k = k)
        }
      } else {# if several candidate variables
        # built the table
        aic.kernel.remove <- data.frame(Part = rep("kernel", length(kernel.remove)), AIC = rep(NA, length(kernel.remove)))
        # rownames of the table are candidate variables
        row.names(aic.kernel.remove) <- paste0("- ", kernel.remove)
        # compute model and AIC for each candidate variable
        for (i in 1:length(kernel.remove)) {
          # new kernel formula
          formula.new <- update.formula(old = m0$kernel.info$Ker1$kernel.formula, new = as.formula(paste0("~ . - ", kernel.remove[i])))
          # new model
          fn <- as.formula(eval(paste0("~Kernel(",paste(formula.new, collapse = ""),", scale = ", kernel.scale, ", kernel.function = '", kernel.function, "', rho = ", kernel.rho, ", gamma = ", kernel.gamma, ", d = ", kernel.d, ")")))
          m.new <- kspm(response = object$Y, linear = m0$linear.formula, kernel = fn , data = dfrm, level = 0, control = kernel.control)
          model.remove.kernel[[i]] <- m.new
          # AIC
          aic.kernel.remove[i, "AIC"] <- extractAIC.kspm(m.new, k = k)
        }
      }

      #-------------------------------------------------------------------------------
      # SUMMARY OF POSSIBLE MODELS
      #-------------------------------------------------------------------------------

      # table of current and new candidate models with AIC
      aic.current <- data.frame(Part = "", AIC = aic)
      row.names(aic.current) <- "<none>"
      aic.new <- rbind(aic.linear.add, aic.kernel.add, aic.linear.remove, aic.kernel.remove)
      aic.new.and.current <- rbind(aic.current, aic.linear.add, aic.kernel.add, aic.linear.remove, aic.kernel.remove)
      aic.new.and.current.order <- aic.new.and.current[order(aic.new.and.current$AIC), ]


      # print information (if asked for)
      if (trace) {
        if (step > 0) {
          cat(paste0("\n\nStep: AIC = ", round(aic,1)))
          cat(paste0("\nLinear part: ~ ", m0$linear.formula)[2])
          cat(paste0("\nKernel part: ~ ", m0$kernel.info$Ker1$kernel.formula)[2])
          cat("\n\n")
        }
        print(aic.new.and.current.order)
      }


      # add information to choose the best model
      aic.new$add.or.remove <- ifelse(substr(row.names(aic.new), 1, 1) == "+", "add", "remove")
      aic.min <- min(aic.new$AIC, na.rm = TRUE)
      aic.new.min <- aic.new[aic.new$AIC == aic.min, ]

      #-------------------------------------------------------------------------------
      # DEFINE THE BEST MODEL
      #-------------------------------------------------------------------------------

      ### if no better or equivalent model, stop the procedure (aic < best.aic.new)
      if (aic < aic.min) {
        continue <- FALSE
      }

      ### if no better model but an equivalent model (aic = best.aic.new)
      if (aic == aic.min) {
        # if only one possibility, continue only if consists in removing a variable
        aic.new.min.remove <- aic.new[aic.new$AIC == aic.min & aic.new$add.or.remove == "remove", ]
        if (dim(aic.new.min.remove)[1] == 0) {
          continue <- FALSE
        } else {# if remove, remove preferably a candidate variable in linear part
          part.to.change <- sort(aic.new.min.remove$Part, decreasing = TRUE)[1]
          if (part.to.change == "linear") {
            i <- which(aic.linear.remove$AIC == aic.min)[1]
            m0 <- model.remove.linear[[i]]
          } else {# if kernel
            i <- which(aic.kernel.remove$AIC == aic.min)[1]
            m0 <- model.remove.kernel[[i]]
          }
        }
      }

      ### if a better model exist (aic > best.aic.new)
      if (aic > aic.min) {
        # if only one solution
        if (dim(aic.new.min)[1] == 1) {
          part.to.change <- aic.new.min$Part[1]
          type.of.change <- aic.new.min$add.or.remove[1]
        } else {
          # if several solutions, the preference order is 1-linear/remove, 2-kernel/remove, 3-kernel/add, 4-linear/add.
          aic.new.min$preference <- ifelse(aic.new.min$Part == "linear" & aic.new.min$add.to.remove == "remove", 1, ifelse(aic.new.min$Part == "kernel" & aic.new.min$add.to.remove == "remove", 2, ifelse(aic.new.min$Part == "kernel" & aic.new.min$add.to.remove == "add", 3, ifelse(aic.new.min$Part == "linear" & aic.new.min$add.to.remove == "add", 4, NA))))
          aic.new.min.ordered <- aic.new.min[order(aic.new.min$preference), ]
          part.to.change <- aic.new.min.ordered$Part[1]
          type.of.change <- aic.new.min.ordered$add.or.remove[1]
        }
        # select the model according to the parameters of the best solution
        if (part.to.change == "linear") {
          if (type.of.change == "remove") {
            i <- which(aic.linear.remove$AIC == aic.min)[1]
            m0 <- model.remove.linear[[i]]
          } else {
            i <- which(aic.linear.add$AIC == aic.min)[1]
            m0 <- model.add.linear[[i]]
          }
        } else {
          if (type.of.change == "remove") {
            i <- which(aic.kernel.remove$AIC == aic.min)[1]
            m0 <- model.remove.kernel[[i]]
          } else {
            i <- which(aic.kernel.add$AIC == aic.min)[1]
            m0 <- model.add.kernel[[i]]
          }
        }
      }
    }
    step <- step + 1
    aic <- extractAIC.kspm(m0, k = k)
  }

  # adapt the call component of the out kspm object
  fn.m0 <- as.formula(eval(paste0("~Kernel(",paste(m0$kernel.info$Ker1$kernel.formula, collapse = ""),", scale = ", kernel.scale, ", kernel.function = '", kernel.function, "', rho = ", kernel.rho, ", gamma = ", kernel.gamma, ", d = ", kernel.d, ")")))

  call.new <- paste0("kspm(response = ", object$call[call.order[3]], ", linear = ", m0$linear.formula[1], m0$linear.formula[2], ", kernel = ", fn.m0[1], fn.m0[2], ", data = ", call.step[call.order.step[5]], ")")
  call.new <- gsub("\"", "'", call.new)
  call.new <- gsub("= ,", "= NULL,", call.new)
  call.new <- gsub("= )", "= NULL)", call.new)
  m0$call <- call.new
  return(invisible(m0))
}


