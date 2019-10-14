#' @title Predicting Kernel Semi parametric Model Fits
#'
#' @description predict method for class "kspm".
#'
#'
#' @param object an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param newdata.linear should be a data frame or design matrix of variables used in the linear part
#' @param newdata.kernel a list containing data frame or design matrix of variables used in each kernel part depending on the specification format of each kernel. When a kernel has been specified using \code{kernel.function = "gram.matrix"} in \code{Kernel} function, the user should also provide the Gram matrix associated to the new data points in \code{newdata.kernel}. The function \link{info.kspm} may help to correctly specify it.
#' @param interval type of interval calculation. If \code{"none"} (default), no interval is computed, if \code{"confidence"}, the confidence interval is computed, if \code{"prediction"}, the prediction interval is computed.
#' @param level confidence level. Default is \code{level = 0.95} meaning 95\% confidence/prediction interval.
#' @param ... further arguments passed to or from other methods.
#'
#' @details \code{predict.kspm} produces predicted values. If a new dataset is not specified, it will return the fitted values from the original data (complete data used in the model specification). If \code{predict.kspm} is applied to a new dataset, all variables used in the original model should be provided in  \code{newdata.linear} and \code{newdata.kernel} arguments but only complete data may be provided. Setting \code{interval} specifies computation of confidence or prediction intervals at the specified \code{level}.
#'
#'
#' @return \code{predict.kspm} returns a vector of predictions or a matrix containing the following components if \code{interval} is set:
#'   \item{fit}{predictions.}
#'   \item{lwr}{lower bound of confidence/prediction intervals.}
#'   \item{upr}{upper bound of confidence/prediction intervals.}
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @importFrom stats model.matrix qt predict
#'
#' @seealso \link{kspm}, \link{summary.kspm}.
#'
#' @examples
#' x <- 1:15
#' z1 <- runif(15, 1, 6)
#' z2 <- rnorm(15, 1, 2)
#' y <- 3*x + (z1 + z2)^2 + rnorm(15, 0, 2)
#' fit <- kspm(y, linear = ~ x, kernel = ~ Kernel(~ z1 + z2,
#' kernel.function = "polynomial", d= 2, rho = 1, gamma = 0))
#' predict(fit, interval = "confidence")
#'
#' @rdname predict.kspm
#' @export predict.kspm
#' @export



predict.kspm <- function(object, newdata.linear = NULL, newdata.kernel = NULL, interval = "none", level = 0.95, ...)
{

  ################################################################################
  # USEFUL INFORMATIONS
  ################################################################################

  # number of samples used in the model
  n <- object$n


  ################################################################################
  # COMPUTE PREDICTION ON INITIAL DATA
  ################################################################################

  # if newdata is null, prediction will be done on data used in the model
  if (is.null(newdata.linear) & is.null(newdata.kernel))
  {
    # fitted values
    fit <- as.matrix(object$fitted.values, ncol = 1)
    # hat matrix for these data
    H <- object$Hat
    # names
    sample.names <- names(object$fitted)

  } else {# if newdata is non null, prediction will be done on this newdata (need to compute new X and new K)

    ################################################################################
    # COMPUTE PREDICTION ON A NEW DATA SET
    ################################################################################

    # new X
    if (is.null(newdata.linear) & dim(object$X)[2] == 1) {
      if (colnames(object$X) == "(Intercept)") {
        X <- matrix(rep(1, dim(newdata.kernel[[1]])[1]), ncol = 1)
      } else {
        stop("missing value for newdata.linear")
      }
    } else {
      X <- model.matrix(object$linear.formula, as.data.frame(newdata.linear))
    }
    # new K
    K <- list()
    # apply kernel function for each main kernel
    for (k in names(object$kernel.info)) {
      # it may be computed if kernel.function is not user.matrix
      if (object$kernel.info[[k]]$type.object != "Gram matrix") {

        if (object$kernel.info[[k]]$type.object == "formula") {
          newdata.Z <- model.frame(object$kernel.info[[k]]$kernel.formula, as.data.frame(newdata.kernel[[k]]))
        }
        if (object$kernel.info[[k]]$type.object == "vector") {
          newdata.Z <- as.matrix(newdata.kernel[[k]], ncol = 1)
        }
        if (object$kernel.info[[k]]$type.object == "design matrix") {
          newdata.Z <- newdata.kernel[[k]]
        }
        # newdata.Z is scaled if necessary
        if (object$kernel.info[[k]]$kernel.scale) {
          newdata.Z.scale <-  scale(newdata.Z, center = object$kernel.info[[k]]$kernel.mean, scale = object$kernel.info[[k]]$kernel.sd)
          Z <- rbind(object$kernel.info[[k]]$Z.scale, newdata.Z.scale)
        } else {
          Z <- rbind(object$kernel.info[[k]]$Z, newdata.Z)
        }
        merge.K <- kernel.matrix(Z = Z, whichkernel = object$kernel.info[[k]]$kernel.function, rho = object$kernel.info[[k]]$rho, gamma = object$kernel.info[[k]]$gamma, d = object$kernel.info[[k]]$d)
        K[[k]] <- merge.K[(n + 1):(dim(merge.K)[2]), 1:n]
      } else {# if the kernel has been defined using user.matrix, the new K is newdata.kernel
        K[[k]] <- newdata.kernel[[k]]
      }
    }
    # compute interaction kernel matrix
    for (k in names(object$lambda)) {
      if (!(k %in% names(object$kernel.info))) {
        inter.k <- strsplit(k, ":")[[1]]
        Km <- K[[inter.k[1]]]
        for (i in 2:length(inter.k)) {
          Km <- Km * K[[inter.k[i]]]
        }
        K[[k]] <- Km
      }
    }

    # fitted values
    fitted.Kalpha <- list()
    KnewGinvM <- list()
    for (k in names(object$lambda)) {
      fitted.Kalpha[[k]] <- K[[k]] %*% object$kernel.coefficients[, k]
      KnewGinvM[[k]] <- K[[k]] %*% object$GinvM[[k]]
    }
    sum_Kalpha <- Reduce("+", fitted.Kalpha)
    if (is.null(sum_Kalpha)) {
      fit <- X %*% object$linear.coefficients
    } else {
      fit <- X %*% object$linear.coefficients + sum_Kalpha
    }

    # hat matrix for a new point
    sum_KGinvM <- Reduce("+", KnewGinvM)
    if (is.null(sum_KGinvM)) {
      H <- X %*% object$XLX_inv %*% t(object$X) %*% object$L
    } else {
      H <- (X - sum_KGinvM %*% object$X) %*% object$XLX_inv %*%  t(object$X) %*% object$L + sum_KGinvM
    }

    # names
    if (!is.null(newdata.linear)) {
      sample.names <- rownames(newdata.linear)
    } else if (!is.null(newdata.kernel)) {
      sample.names <- rownames(newdata.kernel[[1]])
    }  else {
      sample.names <- paste0(1:(dim(X)[1]))
    }

  }


  ################################################################################
  # RETURN FIT INTO VECTOR OR DATA.FRAME
  ################################################################################


  # return vector if no interval
  if (interval == "none") {
    out <- fit[,1]
    names(out) <- sample.names

  } else{# return data.frame

    ################################################################################
    # COMPUTE CONFIDENCE/PREDICTION INTERVALS
    ################################################################################

    # variance-covariance matrix
    fitted.vcov <- H %*% (object$sigma * diag(1, ncol = n, nrow = n)) %*% t(H)

    # variance of fitted values
    fitted.var <- diag(fitted.vcov)
    # SE of fitted values
    fitted.se <- sqrt(fitted.var)
    quantile <- -qt((1 - level)/2, object$edf)
    # confidence interval
    if (interval == "confidence") {
      lwr <- fit - quantile*fitted.se
      upr <- fit + quantile*fitted.se
      out <- data.frame(fit = fit, lwr = lwr, upr = upr)
      rownames(out) <- sample.names
    }
    # prediction interval
    if (interval == "prediction") {
      lwr <- fit - quantile*sqrt(fitted.se^2 + object$sigma^2)
      upr <- fit + quantile*sqrt(fitted.se^2 + object$sigma^2)
      out <- data.frame(fit = fit, lwr = lwr, upr = upr)
      rownames(out) <- sample.names
    }
  }

  # return results
  return(out)
}
