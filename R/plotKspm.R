#' @title Plot Diagnostics for a kspm Object
#'
#' @description Five plots (selectable by \code{which}) are currently available: a plot of residuals against fitted values, a scale Location plot of \eqn{\sqrt{\mid residuals \mid}}{sqrt(| residuals |)} against fitted values, a Normal Q Q plot for residuals, a plot of Cook's distances versus row labels and a plot of residuals against leverages. By default, the first three and 5 are provided.
#'
#'
#' @param x an object of class "kspm", usually, a result of a call to \code{kspm}.
#' @param which if a subset of the plots is required, specify a subset of the numbers 1:5.
#' @param cook.levels levels of Cook's distance at which to draw contours.
#' @param id.n number of points to be labelled in each plot, starting with the most extreme.
#' @param labels.id vector of labels, from which the labels for extreme points will be chosen. NULL uses names associated to response specified in \code{kspm}.
#' @param cex.id size of point labels.
#' @param col.id color of point labels.
#' @param ... further arguments passed to or from other methods.
#'
#'
#' @author Catherine Schramm, Aurelie Labbe, Celia Greenwood
#'
#' @references
#'
#' Kim, Choongrak, Byeong U. Park, and Woochul Kim. "Influence diagnostics in semiparametric regression models." Statistics and probability letters 60.1 (2002): 49:58.
#'
#' @examples
#' x <- 1:15
#' z1 <- runif(15, 1, 6)
#' z2 <- rnorm(15, 1, 2)
#' y <- 3*x + (z1 + z2)^2 + rnorm(15, 0, 2)
#' fit <- kspm(y, linear = ~ x, kernel = ~ Kernel(~ z1 + z2,
#' kernel.function = "polynomial", d= 2, rho = 1, gamma = 0))
#' plot(fit)
#'
#' @seealso \link{kspm} for fitting the model, \link{summary.kspm} for summary
#'
#' @importFrom graphics segments text
#' @importFrom stats qqline qqnorm
#'
#' @export
#'
#'






plot.kspm <- function(x, which = c(1:3, 5), cook.levels = c(0.5, 1.0), id.n = 3, labels.id = names(x$residuals), cex.id = 0.75, col.id = "blue", ...)
{

  ################################################################################
  # COMPUTE STANDARDIZED RESIDUALS
  ################################################################################

  leverage <- diag(x$Hat)
  standardized.residuals <- x$residuals / (x$sigma * sqrt(1 - leverage))

  ################################################################################
  # RESIDUALS AGAINST FITTED VALUES (which = 1)
  ################################################################################

  if (1 %in% which) {
    a <- max(abs(x$residuals))
    plot(x$fitted, x$residuals, xlab = "Fitted", ylab = "Residuals", main = "Residuals vs Fitted", ylim = c(-a, a), ...)
    abline(h = 0, col = "gray", lty = 2)
    # most extreme values (absolute value of residual)
    id.labels <- names(sort(abs(x$residuals), decreasing = TRUE))[1:id.n]
    for (id in id.labels) {
      text(x$fitted[id], x$residuals[id], id, pos = 3, col = col.id, cex = cex.id, xpd = TRUE)
    }
    # ask for new graphic
    devAskNewPage(ask = TRUE)
  }

  ################################################################################
  # SCALE-LOCATION AGAINST FITTED VALUES (which = 2)
  ################################################################################

  if (2 %in% which) {
    plot(x$fitted, sqrt(abs(standardized.residuals)), xlab = "Fitted", ylab = expression(sqrt(abs("Standardized residuals"))), main = "Scale-Location", ...)
    # most extreme values
    id.labels <- names(sort(abs(standardized.residuals), decreasing = TRUE))[1:id.n]
    for (id in id.labels) {
      text(x$fitted[id], sqrt(abs(standardized.residuals))[id], id, pos = 3, col = col.id, cex = cex.id, xpd = TRUE)
    }
    # ask for new graphic
    devAskNewPage(ask = TRUE)
  }

  ################################################################################
  # NORMAL Q-Q PLOT (which = 3)
  ################################################################################

  if (3 %in% which) {
    qq <- qqnorm(standardized.residuals, ylab = "Standardized residuals", ...)
    qqline(standardized.residuals, col = "gray", lty = 2)
    # most extreme values
    id.labels <- names(sort(abs(standardized.residuals), decreasing = TRUE))[1:id.n]
    for (id in id.labels) {
      text(qq$x[which(names(standardized.residuals) == id)], qq$y[which(names(standardized.residuals) == id)], id, pos = 3, col = col.id, cex = cex.id, xpd = TRUE)
    }
    # ask for new graphic
    devAskNewPage(ask = TRUE)
  }

  ################################################################################
  # COOK'S DISTANCE (which = 4)
  ################################################################################

  if (4 %in% which) {
    cook.distance <- x$residuals ^ 2 * leverage / (x$sigma ^ 2 * sum(leverage) * (1 - leverage) ^ 2)
    names(cook.distance) <- labels.id
    plot(0, 0, col = "white", main = "Cook's distance", xlim = c(0, x$n), ylim = c(0, max(cook.distance)), xlab = "Obs. number", ylab = "Cook's distance", ...)
    segments(1:x$n, 0, 1:x$n, cook.distance)
    # most extreme values
    id.labels <- names(sort(cook.distance, decreasing = TRUE))[1:id.n]
    for (id in id.labels) {
      text(which(names(cook.distance) == id), cook.distance[id], id, pos = 3, col = col.id, cex = cex.id, xpd = TRUE)
    }
    # ask for new graphic
    devAskNewPage(ask = TRUE)
  }

  ################################################################################
  # STANDARDIZED RESIDUALS AGAINST LEVERAGE (which = 5)
  ################################################################################

  if (5 %in% which) {
    plot(leverage, standardized.residuals, xlim = c(0, max(leverage)), main = "Residuals vs Leverage", xlab = "Leverage", ylab = "Standardized residuals", ...)
    names(leverage) <- labels.id
    abline(h = 0, col = "gray", lty = 3)
    abline(v = 0, col = "gray", lty = 3)
    # most extreme values
    cook.distance <- x$residuals ^ 2 * leverage / (x$sigma ^ 2 * sum(leverage) * (1 - leverage) ^ 2)
    names(cook.distance) <- labels.id
    id.labels <- names(sort(abs(cook.distance), decreasing = TRUE))[1:id.n]
    for (id in id.labels) {
      text(leverage[id], standardized.residuals[id], id, pos = 3, col = col.id, cex = cex.id, xpd = TRUE)
    }
    for (i in cook.levels) {
      leverage.sort <- sort(leverage)
      cook.pos <- sqrt(i * sum(leverage.sort) * (1 - leverage.sort) / leverage.sort)
      lines(leverage.sort, cook.pos, col = "red", lty = 2)
      text(leverage.sort[x$n],cook.pos[x$n], i, pos = 4, col = "red", xpd = TRUE)
      cook.neg <- -sqrt(i * sum(leverage.sort) * (1 - leverage.sort) / leverage.sort)
      lines(leverage.sort, cook.neg, col = "red", lty = 2)
      text(leverage.sort[x$n],cook.neg[x$n], i, pos = 4, col = "red", xpd = TRUE)
    }
    legend("bottomleft", legend = "Cook's distance", col = "red", lty = 2, bty = "n")
    # ask for new graphic
    devAskNewPage(ask = TRUE)
  }

  #-------------------------------------------------------------------------------
  # INFLUENCE FOR LINEAR PART
  #-------------------------------------------------------------------------------

  if (sum(c(6, 7, 8, 9) %in% which) > 0) {
    H.linear <- x$X %*% x$XLX_inv %*% t(x$X) %*% x$L
    leverage.linear <- diag(H.linear)
    residuals.linear <- ((diag(1, nrow = x$n, ncol = x$n) - H.linear) %*% x$Y)[, 1]
    standardized.residuals.linear <- residuals.linear / (x$sigma * sqrt(1 - leverage.linear))
  }

  ################################################################################
  # LINEAR PART: COOK'S DISTANCE (which = 6)
  ################################################################################

  if (6 %in% which) {
    cook.distance.linear <- residuals.linear ^ 2 * leverage.linear / (x$sigma ^ 2 * sum(leverage.linear) * (1 - leverage.linear) ^ 2)
    names(cook.distance.linear) <- labels.id
    plot(0, 0, col = "white", main = "Cook's distance (Linear part)", xlim = c(0, x$n), ylim = c(0, max(cook.distance.linear)), xlab = "Obs. number", ylab = "Cook's distance", ...)
    segments(1:x$n, 0, 1:x$n, cook.distance.linear)
    # most extreme values
    id.labels <- names(sort(cook.distance.linear, decreasing = TRUE))[1:id.n]
    for (id in id.labels) {
      text(which(names(cook.distance.linear) == id), cook.distance.linear[id], id, pos = 3, col = col.id, cex = cex.id, xpd = TRUE)
    }
    # ask for new graphic
    devAskNewPage(ask = TRUE)
  }

  ################################################################################
  # LINEAR PART: STANDARDIZED RESIDUALS AGAINST LEVERAGE (which = 7)
  ################################################################################

  if (7 %in% which) {
    plot(leverage.linear, standardized.residuals.linear, xlim = c(0, max(leverage.linear)), main = "Residuals vs Leverage (Linear part)", xlab = "Leverage", ylab = "Standardized residuals", ...)
    names(leverage.linear) <- labels.id
    abline(h = 0, col = "gray", lty = 3)
    abline(v = 0, col = "gray", lty = 3)
    # most extreme values
    cook.distance.linear <- residuals.linear ^ 2 * leverage.linear / (x$sigma ^ 2 * sum(leverage.linear) * (1 - leverage.linear) ^ 2)
    names(cook.distance.linear) <- labels.id
    id.labels <- names(sort(abs(cook.distance.linear), decreasing = TRUE))[1:id.n]
    for (id in id.labels) {
      text(leverage.linear[id], standardized.residuals.linear[id], id, pos = 3, col = col.id, cex = cex.id, xpd = TRUE)
    }
    for (i in cook.levels) {
      leverage.sort <- sort(leverage.linear)
      cook.pos <- sqrt(i * sum(leverage.sort) * (1 - leverage.sort) / leverage.sort)
      lines(leverage.sort, cook.pos, col = "red", lty = 2)
      text(leverage.sort[x$n],cook.pos[x$n], i, pos = 4, col = "red", xpd = TRUE)
      cook.neg <- -sqrt(i * sum(leverage.sort) * (1 - leverage.sort) / leverage.sort)
      lines(leverage.sort, cook.neg, col = "red", lty = 2)
      text(leverage.sort[x$n],cook.neg[x$n], i, pos = 4, col = "red", xpd = TRUE)
    }
    legend("bottomleft", legend = "Cook's distance", col = "red", lty = 2, bty = "n")
    # ask for new graphic
    devAskNewPage(ask = TRUE)
  }



  #-------------------------------------------------------------------------------
  # INFLUENCE FOR KERNEL PART
  #-------------------------------------------------------------------------------

  if (sum(c(8, 9) %in% which) > 0) {
    H.kernel <- diag(1, nrow = x$n, ncol = x$n) - x$L
    leverage.kernel <- diag(H.kernel)
    residuals.kernel <- ((diag(1, nrow = x$n, ncol = x$n) - H.kernel) %*% x$Y)[, 1]
    standardized.residuals.kernel <- residuals.kernel / (x$sigma * sqrt(1 - leverage.kernel))
  }

  ################################################################################
  # KERNEL PART: COOK'S DISTANCE (which = 8)
  ################################################################################

  if (8 %in% which) {
    cook.distance.kernel <- (residuals.kernel ^ 2) * (leverage.kernel ^ 2) / (x$sigma ^ 2 * sum(leverage.kernel) * (1 - leverage.kernel) ^ 2)
    names(cook.distance.kernel) <- labels.id
    plot(0, 0, col = "white", main = "Cook's distance (Kernel part)", xlim = c(0, x$n), ylim = c(0, max(cook.distance.kernel)), xlab = "Obs. number", ylab = "Cook's distance", ...)
    segments(1:x$n, 0, 1:x$n, cook.distance.kernel)
    # most extreme values
    id.labels <- names(sort(cook.distance.kernel, decreasing = TRUE))[1:id.n]
    for (id in id.labels) {
      text(which(names(cook.distance.kernel) == id), cook.distance.kernel[id], id, pos = 3, col = col.id, cex = cex.id, xpd = TRUE)
    }
    # ask for new graphic
    devAskNewPage(ask = TRUE)
  }

  ################################################################################
  # LINEAR PART: STANDARDIZED RESIDUALS AGAINST LEVERAGE (which = 9)
  ################################################################################

  if (9 %in% which) {
    plot(leverage.kernel, standardized.residuals.kernel, xlim = c(0, max(leverage.kernel)), main = "Residuals vs Leverage (Kernel part)", xlab = "Leverage", ylab = "Standardized residuals", ...)
    names(leverage.kernel) <- labels.id
    abline(h = 0, col = "gray", lty = 3)
    abline(v = 0, col = "gray", lty = 3)
    # most extreme values
    cook.distance.kernel <- (residuals.kernel ^ 2) * (leverage.kernel ^ 2) / (x$sigma ^ 2 * sum(leverage.kernel) * (1 - leverage.kernel) ^ 2)
    names(cook.distance.kernel) <- labels.id
    id.labels <- names(sort(abs(cook.distance.kernel), decreasing = TRUE))[1:id.n]

    for (id in id.labels) {
      text(leverage.kernel[id], standardized.residuals.kernel[id], id, pos = 3, col = col.id, cex = cex.id, xpd = TRUE)
    }
    for (i in cook.levels) {
      leverage.sort <- sort(leverage.kernel)
      cook.pos <- sqrt(i * sum(leverage.sort) * (1 - leverage.sort) / (leverage.sort ^ 2))
      lines(leverage.sort, cook.pos, col = "red", lty = 2)
      text(leverage.sort[x$n],cook.pos[x$n], i, pos = 4, col = "red", xpd = TRUE)
      cook.neg <- -sqrt(i * sum(leverage.sort) * (1 - leverage.sort) / leverage.sort)
      lines(leverage.sort, cook.neg, col = "red", lty = 2)
      text(leverage.sort[x$n],cook.neg[x$n], i, pos = 4, col = "red", xpd = TRUE)
    }
    legend("bottomleft", legend = "Cook's distance", col = "red", lty = 2, bty = "n")
    # ask for new graphic
    devAskNewPage(ask = TRUE)
  }

  # stop asking for new graphic
  devAskNewPage(ask = FALSE)

}



