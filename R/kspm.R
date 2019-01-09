#' @title Fitting Kernel Semi Parametric model
#'
#' @description kspm is used to fit kernel semi parametric models.
#'
#' @param response a character with the name of the response variable or a vector containing the outcome or a matrix with outcome in the first column.
#' @param linear an optional object of class "formula": a symbolic description of the linear part of the model to be fitted or a vector or a matrix containing covariates included in the linear part of the model. Default is intercept only. The details of model specification are given under ‘Details’.
#' @param kernel an object of class "formula": a symbolic description of the kernel part of the model to be fitted. If missing a linear model is fitted using lm function. The details of model specification are given under ‘Details’.
#' @param data  an optional data frame containing the variables in the model. If NULL (default), data are taken from the workspace.
#' @param level printed information about the model (0: no information, 1: information about kernels included in the model (default))
#' @param control see \link{kspmControl}.
#'
#' @details The kernel semi parametric model refers to the following equation \eqn{Y_i = X_i\beta + h(Z_i) + e_i}{Y_i = X_iB + h(Z_i) + e_i} with \eqn{i=1..n}{i=1..n} where \eqn{n}{n} is the sample size, \eqn{Y}{Y} is the univariate response, \eqn{X\beta}{XB} is the linear part, \eqn{h(Z)}{h(Z)} is the kernel part and \eqn{e}{e} are the residuals. The linear part is defined using the \code{linear} argument by specifying the covariates \eqn{X}{X}. It could be either a formula, a vector of length \eqn{n}{n} if only one variable is included in the linear part or a \eqn{n \times p}{n x p} design matrix containing the values of the \eqn{p}{p} covariates included in the linear part (columns), for each individuals (rows). By default, an intercept is included. To remove the intercept term, use formula specification and add the term \code{-1}, as usual. Kernel part is defined using the \code{kernel} argument. It should be a formula of \code{Kernel} object(s). For a multiple kernel semi parametric model, \code{Kernel} objects are separated by the usual signs \code{"+"}, \code{"*"} and \code{":"} to specify addition and interaction between kernels. Specification formats of each \code{Kernel}  object may be different. See \link{Kernel} for more information about their specification.
#'
#'
#' @return \code{kspm} returns an object of class kspm.
#' @return An object of class kspm is a list containing the following components:
#'  \item{linear.coefficients}{matrix of coefficients associated with linear part, the number of coefficients is the number of terms included in linear part}
#'  \item{kernel.coefficients}{matrix of coefficients associated with kernel part, the number of rows is the sample size included in the analysis and the number of columns is the number of kernels included in the model}
#'  \item{lambda}{penalization parameter(s)}
#'  \item{fitted.values}{the fitted mean values}
#'  \item{residuals}{the residuals, that is response minus the fitted values}
#'  \item{sigma}{standard deviation of residuals}
#'  \item{Y}{vector of responses}
#'  \item{X}{design matrix for linear part}
#'  \item{K}{kernel matrices computed by the model}
#'  \item{n.total}{total sample size}
#'  \item{n}{sample size of the model (model is performed on complete data only)}
#'  \item{edf}{effective degree of freedom}
#'  \item{linear.formula}{formula corresponding to the linear part of the model}
#'  \item{kernel.info}{information about kernels included in the model such as matrices of covariates (\code{Z}),  kernel function (\code{type}), values of hyperparameters (\code{rho}, \code{gamma}, \code{d}). A boolean indicates if covariates were scaled (\code{kernel.scale}) and if \code{TRUE}, \code{kernel.mean}, \code{kernel.sd} and \code{Z.scale} give information about scaling. \code{kernel.formula} indicates the formula of the kernel and \code{free.parameters} indicates the hyperparameters that were estimated by the model.}
#'  \item{Hat}{The hat matrix \eqn{H}{H} such that \eqn{\hat{Y} = HY}{Y_hat = HY}}
#'  \item{L}{A matrix corresponding to \eqn{I - \sum\limits_{\ell = 1}^L K_{\ell} G_{\ell}^{-1} M_{\ell}}{I-sum_l K_l G_l^{-1}L_l} according to our notations}
#'  \item{XLX_inv}{A matrix corresponding to \eqn{(XLX)^{-1}}{(XLX)^{-1}}}
#'  \item{GinvM}{A list of matrix, each corresponding to a kernel and equaling \eqn{G_{\ell}^{-1}M_{\ell}}{G_l^{-1}M_l} according to our notations}
#'  \item{control}{List of control parameters}
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @references
#' Liu, D., Lin, X., and Ghosh, D. (2007). Semiparametric regression of multidimensional genetic pathway data: least squares kernel machines and linear mixed models. Biometrics, 63(4), 1079:1088.
#'
#' Kim, Choongrak, Byeong U. Park, and Woochul Kim. "Influence diagnostics in semiparametric regression models." Statistics and probability letters 60.1 (2002): 49:58.
#'
#' Oualkacha, Karim, et al. "Adjusted sequence kernel association test for rare variants controlling for cryptic and family relatedness." Genetic epidemiology 37.4 (2013): 366:376.
#'
#' @examples
#' x <- 1:15
#' z1 <- runif(15, 1, 6)
#' z2 <- rnorm(15, 1, 2)
#' y <- 3*x + (z1 + z2)^2 + rnorm(15, 0, 2)
#' fit <- kspm(y, linear = ~ x, kernel = ~ Kernel(~ z1 + z2,
#' kernel.function = "polynomial", d= 2, rho = 1, gamma = 0))
#' summary(fit)
#'
#' @importFrom stats setNames
#'
#' @seealso \link{summary.kspm} for summary, \link{predict.kspm} for predictions, \link{plot.kspm} for diagnostics
#'
#' @export




kspm <- function(response, linear = NULL, kernel = NULL, data = NULL, level = 1, control = kspmControl())
{

  ################################################################################
  # CALL
  ################################################################################

  call <- match.call(expand.dots = FALSE)
  call.names <- paste0(call)
  call.order <- match(c("response", "linear", "kernel", "data", "level", "control"), names(call), 0L)

  ################################################################################
  # SOME CHECKS
  ################################################################################

  # data should be a data.frame or a matrix
  if (!is.null(data) && (!inherits(data, "data.frame") & !inherits(data, "matrix"))) {
    stop("data should be a data.frame or a matrix")
  }
  # level should be >= 0
  if (!(level %in% c(0, 1, 2))) {
    stop("level should be 0, 1 or 2")
  }

  ################################################################################
  # NAMES
  ################################################################################

  if (!is.null(data)) {
    names <- rownames(data)
  } else if (is.matrix(response) & !is.null(rownames(response)) & length(rownames(response)) == length(unique(rownames(response)))) {
    names <- rownames(response)
  } else if (is.vector(response) & !is.null(names(response)) & length(names(response)) == length(unique(names(response)))) {
    names <- names(response)
  } else {
    names <- paste0(1:length(response))
  }

  ################################################################################
  # RESPONSE Y
  ################################################################################

  # if no response => STOP
  if (is.null(response)) {
    stop("reponse is missing")
  } else if (is.character(response)) { # if response is character => search corresponding column in data
    if (is.null(data)) {
      stop("response is a character but data is missing")
    } else {
      Y <- matrix(data[, response], nrow = dim(data)[1])
      colnames(Y) <- response
      rownames(Y) <- names
    }
  } else if (is.vector(response)) { # if response is a vector
    if (!is.numeric(response)) {
      stop("response vector should be numeric")
    } else {
      Y <- matrix(response, nrow = length(response))
      colnames(Y) <- call.names[call.order[1]]
      rownames(Y) <- names
    }
  } else if (is.matrix(response)) { # if response is a matrix => consider first column only
    if (dim(response)[2] > 1) {
      warning("several entries for response, only the first column is used")
    }
    Y <- matrix(response[, 1], nrow = dim(response)[1])
    colnames(Y) <- colnames(response)[1]
    rownames(Y) <- names
  } else {
    stop("reponse should be a character, a vector or a matrix")
  }

  ################################################################################
  # LINEAR X
  ################################################################################

  # if linear is NULL => Intercept only
  if (is.null(linear)) {
    X <- matrix(rep(1, dim(Y)[1]), dim(Y)[1])
    colnames(X) <- "(Intercept)"
    rownames(X) <- names
    linear.formula <- ~1
  } else if (inherits(linear, "formula")) { # if linear is a formula
    linear.formula <- linear
    options(na.action = "na.pass")
    if (is.null(data)) {
      X <- model.matrix(linear)
      rownames(X) <- names
    } else {
      X <- model.matrix(linear, data = data)
      rownames(X) <- names
    }
    options(na.action = "na.omit")
  } else if (is.vector(linear)) { # if linear is a vector
    X <- matrix(cbind(rep(1, length(linear)), linear), nrow = length(linear))
    colnames(X) <- c("(Intercept)", call.names[call.order[2]])
    rownames(X) <- names
    linear.formula <- as.formula(paste0("~", call.names[call.order[2]]))
  } else if (is.matrix(linear)) {  # if linear is a matrix
    if (is.null(colnames(linear))) {
      colnames(linear) <- paste0("Lin", 1:(dim(linear)[2]))
    }
    X <- cbind(matrix(rep(1, dim(linear)[1]), dim(linear)[1]), linear)
    colnames(X) <- c("(Intercept)", colnames(linear))
    rownames(X) <- names
    linear.formula <- as.formula(paste0("~", paste(colnames(linear), collapse = " + ")))
  } else {
    stop("linear should be a formula, a vector or a matrix")
  }

  # Combine data from Y and X
  cbind.data <- cbind(Y, X)

  ################################################################################
  # KSPM
  ################################################################################

  if (!is.null(kernel) ) {

    ################################################################################
    # EVALUTE DATA FOR KERNEL Z/K
    ################################################################################
    if (!inherits(kernel, "formula")) {
      stop("kernel should be a formula")
    }

    # information on data used in each kernel
    kernelObject <- objects.Kernel(kernel)[[1]]
    informUserMatrix <- objects.Kernel(kernel)[[2]]
    # for each kernel, evalute the data
    for (i in 1:length(kernelObject)) {
      # data used in the kernel
      obj = eval(parse(text = kernelObject[[i]]), envir = parent.frame())
      # if obj is a formula
      if (inherits(obj, "formula")) {
        options(na.action = "na.pass")
        if (is.null(data)) {
          Z.obj <- model.matrix(obj)[, -1]
        } else {
          Z.obj <- model.matrix(obj, data = data)[, -1]
        }
        options(na.action = "na.omit")
      }
      # if obj is a vector
      if (is.vector(obj)) {
        Z.obj <- matrix(obj, nrow = length(obj))
      }
      # if obj is a matrix
      if (is.matrix(obj)) {

        if (informUserMatrix[i] != "gram.matrix") {
          Z.obj <- obj
        } else {
          K.obj <- obj
          Z.obj <- matrix(ifelse(apply(K.obj, 1, function(x) sum(is.na(x))) == dim(K.obj)[1], NA, 1), nrow = dim(K.obj)[1])
          # compute.kernel <- FALSE
        }
      }
      # add data information to cbind.data
      cbind.data <- cbind(cbind.data, Z.obj)
    }

    ################################################################################
    # REMOVE MISSING DATA
    ################################################################################

    # sample size (including na)
    n.total <- dim(cbind.data)[1]
    # keep not missing data
    not.missing <- rownames(as.data.frame(na.omit(cbind.data)))
    # sample size (excluding na)
    n <- length(not.missing)

    # Y without NA
    Y <- matrix(Y[not.missing, ], nrow = length(not.missing))
    rownames(Y) <- not.missing

    # X without NA
    X.names <- colnames(X)
    X <- matrix(X[not.missing, ], nrow = length(not.missing))
    rownames(X) <- not.missing
    colnames(X) <- X.names

    ################################################################################
    # COMPUTE KERNEL Z/K ON NON MISSING DATA
    ################################################################################

    # kernel data # VERIFY NAMES AND NOT.MISSING
    kernelList <- kernel.list(kernel, data = data, names = names)

    ################################################################################
    # INFORMATION ON WHICH KERNEL SHOULD BE PUT IN THE MODEL
    ################################################################################

    # kernels included in the model
    inclKernel <- kernelList$inclKernel
    # names of kernels
    nameKernel <- kernelList$nameKernel
    singKernel <- kernelList$singKernel
    # list of kernel expression in the current function
    mergExpres <- kernelList$mergExpres

    # print information about kernels included in the model
    if (level > 0) {
      cat("\n\n------------------------------------------------------------------ ")
      cat("\nThe model includes the following kernels: ")
      cat(paste0("\n \n", nameKernel[inclKernel]))
      cat("\n\n------------------------------------------------------------------ ")
      cat("\nDetails: ")
      for (i in 1:length(singKernel)) {
        cat(paste0("\n \n", nameKernel[i], ": ", mergExpres[i]))
        if (is.list(mergExpres[[i]])) {
          cat(paste0("\n", nameKernel[i], " has been expressed with ", length(mergExpres[[i]]), " different expressions representing the same kernel"))
        }
      }
      cat("\n\n------------------------------------------------------------------\n ")
    }

    ################################################################################
    # COMPUTE MODEL
    ################################################################################

    # get coefficients and lambda parameter
    param <- search.parameters(Y = Y, X = X, kernelList = kernelList, n = n, not.missing = not.missing, compute.kernel = TRUE, controlKspm = control)

    # fitted values
    fitted.Kalpha <- list()
    for (i in 1:length(param$K)) {

      fitted.Kalpha[[i]] <- param$K[[i]] %*% param$alpha[[i]]
    }
    if (dim(X)[2] == 0) {
      fitted.values <- Reduce("+", fitted.Kalpha)
    } else {
      fitted.values <- X %*% param$beta + Reduce("+", fitted.Kalpha)
    }

    # residuals
    residuals <- Y - fitted.values

    # standard deviation of residuals
    edf = n - sum(diag(param$Hat))
    sigma <- sqrt((1 / edf) * sum(residuals ^ 2))

    ################################################################################
    # DEFINITION KSPM
    ################################################################################

    # Informations about kernels
    kernelList.single = list()
    for (i in 1:length(param$kernelList$singKernel)) {
      kernelList.single[[i]] = param$kernelList$listKernel[[i]][c("Z", "kernel.function", "rho", "gamma", "d", "kernel.scale", "kernel.mean", "kernel.sd", "Z.scale", "kernel.formula", "free.parameters", "type.object")]
    }
    names(kernelList.single) = kernelList$nameKernel[kernelList$singKernel]

    out <- list(

      # Call
      call = call,

      # Linear coefficients
      linear.coefficients = param$beta,
      # Kernel coefficients
      kernel.coefficients = matrix(unlist(param$alpha), nrow = n, dimnames = list(not.missing, nameKernel[inclKernel])),

      # Lambda
      lambda = setNames(param$lambda, nameKernel[inclKernel]),

      # Fitted values
      fitted.values = fitted.values[, 1],

      # Residuals
      residuals = residuals[, 1],
      # Standard deviation of residuals
      sigma = sigma,

      # Response matrix
      Y = Y[, 1],

      # Linear part matrix
      X = X,

      # Kernel matrix (computed on scaled Z matrix if scale was applied)
      K = setNames(param$K, nameKernel[inclKernel]),

      # Hat matrix
      Hat = param$Hat,
      # L matrix
      L = param$L,
      # XLX_inv matrix
      XLX_inv = param$XLX_inv,
      # list of GinvM matrix
      GinvM = setNames(param$GinvM, nameKernel[inclKernel]),

      # Original sample size
      n.total = n.total,
      # Sample size used by the model (NA excluded)
      n = n,
      # Effective degrees of freedom of the model (n - tr(H))
      edf = edf,

      # Formula for linear part
      linear.formula = linear.formula,
      # Informations about kernel(s)
      kernel.info = kernelList.single,
      # Information about control
      control = control)
  }

  ################################################################################
  # LM
  ################################################################################

  if (is.null(kernel)) {

    if (level > 0) {
      cat("No value for kernel, a lm model is performed instead")
    }


    ################################################################################
    # CHECK DIMENSIONS Y, X and ROWNAMES
    ################################################################################

    # Check compatibility of dimensions between Y, X
    if (dim(Y)[1] != dim(X)[1]) {
      stop("the lengths of the variables differ")
    }

    # Check compatibility of rownames between Y, X
    if ( !is.null(rownames(Y)) && !is.null(rownames(X)) ) {
      if (sum(rownames(Y) != rownames(X)) > 0) {
        stop("sample names differ or are not in the same order")
      }
    }

    ################################################################################
    # REMOVE MISSING DATA Y, X
    ################################################################################

    cbind.data <- cbind(Y, X)
    if (is.null(rownames(cbind.data))) {
      rownames(cbind.data) <- 1:dim(Y)[1]
    }

    cbind.data.excludeNA <- as.data.frame(na.omit(cbind.data))
    if (dim(cbind.data.excludeNA)[1] != dim(cbind.data)[1]) {
      # Y without NA
      Y <- matrix(cbind.data.excludeNA[,1], ncol = 1)
      colnames(Y) <- colnames(cbind.data.excludeNA)[1]
      rownames(Y) <- rownames(cbind.data.excludeNA)
      # X without NA
      X <- as.matrix(cbind.data.excludeNA[,2:dim(cbind.data)[2]], ncol = dim(cbind.data)[2] - 1)
      colnames(X) <- colnames(cbind.data.excludeNA)[2:(dim(cbind.data)[2])]
      rownames(X) <- rownames(cbind.data.excludeNA)
    }

    ################################################################################
    # SAMPLE SIZE
    ################################################################################

    # sample size (including na)
    n.total <- dim(cbind.data)[1]
    # sample size (excluding na)
    n <- dim(cbind.data.excludeNA)[1]

    ################################################################################
    # COMPUTE MODEL
    ################################################################################

    # fit lm
    fit.lm <- lm(Y ~ X - 1)

    # coefficients with correct names
    lin.coef <- as.matrix(fit.lm$coefficients)
    new.names <- substr(rownames(as.matrix(lin.coef)), 2, nchar(rownames(as.matrix(lin.coef))))
    row.names(lin.coef) <- new.names


    # hat matrix (H)
    Hat <- X %*% solve(t(X) %*% X) %*% t(X)

    ################################################################################
    # DEFINITION LM OUTPUT
    ################################################################################

    # ATTENTION MODIFIER FORMULA OUTPUT
    out <- list(

      # Call
      call = call,

      # Linear coefficients
      linear.coefficients = lin.coef,
      # Kernel coefficients
      kernel.coefficients = NULL,

      # Lambda
      lambda = NULL,

      # Fitted values
      fitted.values = fit.lm$fitted.values,

      # Residuals
      residuals = residuals(fit.lm),
      # Standard deviation of residuals
      sigma = summary(fit.lm)$sigma,

      # Response matrix
      Y = Y[, 1],

      # Linear part matrix
      X = X,

      # Kernel matrix (computed on scaled Z matrix if scale was applied)
      K = NULL,

      # Hat matrix
      Hat = Hat,
      # L matrix
      L = diag(1, n),
      # XLX_inv matrix
      XLX_inv = solve(t(X) %*% X),
      # list of GinvM matrix
      GinvM = diag(1, n), # A VERIFIER

      # Original sample size
      n.total = n.total,
      # Sample size used by the model (NA excluded)
      n = n,
      # Effective degrees of freedom of the model (n - tr(H))
      edf = summary(fit.lm)$df[2],

      # Formula for linear part
      linear.formula = linear.formula,
      # Informations about kernel(s)
      kernel.info = NULL,
      # Information about control
      control = control)
  }

  ################################################################################
  # RETURN OUTPUT
  ################################################################################

  # return kspm object
  class(out) <- c("kspm")
  return(out)

}
