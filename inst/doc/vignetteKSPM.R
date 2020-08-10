## ----setup, include = FALSE----------------------------------------------
# Sys.setenv(hydat_eval = "local")
variablelocal <- Sys.getenv("hydat_eval")
LOCAL <- variablelocal == "local"
print(LOCAL)

## ---- eval=LOCAL---------------------------------------------------------
#  
#  
#  library(KSPM)
#  
#  

## ---- eval=LOCAL---------------------------------------------------------
#  data(csm)
#  head(csm)
#  

## ---- eval=LOCAL---------------------------------------------------------
#  csm.fit1 <- kspm(response = "Ratings", kernel = ~Kernel(~Gross + Budget + Screens + Sequel, kernel.function = "gaussian", rho = 61.22), data = csm)

## ---- eval=LOCAL---------------------------------------------------------
#  summary(csm.fit1)

## ---- fig.height=8, fig.show='hold', fig.width=8, eval=LOCAL-------------
#  par(mfrow = c(2,2), mar = c(5, 5, 5, 2))
#  plot(csm.fit1, which = c(1, 3, 5), cex.lab = 1.5, cex = 1.3)
#  hist(csm$Ratings, cex.lab = 1.5, main = "Histogram of Ratings", xlab = "Ratings")

## ---- fig.height=4, fig.show='hold', fig.width=8, eval=LOCAL-------------
#  par(mfrow = c(1,2), mar = c(5,5,5,2))
#  plot(derivatives(csm.fit1), subset = "Gross", main = "Pointwise derivatives according to Gross Income", cex.main = 0.8)
#  plot(derivatives(csm.fit1), subset = "Screens", col = csm$Sequel, pch = 16, main = "Pointwise derivatives according to \n Number of Screens and Sequel", cex.main = 0.8, ylim = c(-0.4, 0.8))
#  legend("topleft", fill = palette()[1:7], legend = 1:7, title = "Sequel", horiz = TRUE, cex = 0.7)

## ---- eval=LOCAL---------------------------------------------------------
#  csm.fit2 <- kspm(response = "Ratings", kernel = ~Kernel(~Gross + Budget + Screens + Sequel, kernel.function = "polynomial", rho = 1, gamma = 1, d = 2), data = csm, level = 0)

## ---- eval=LOCAL---------------------------------------------------------
#  extractAIC(csm.fit1)
#  extractAIC(csm.fit2)

## ---- eval = FALSE-------------------------------------------------------
#  
#  csm.fit3 <- kspm(response = "Ratings", linear = NULL, kernel = ~Kernel(~ Gross + Budget + Screens + Sequel, kernel.function = "gaussian", rho = 61.22) * Kernel(~ Sentiment + Views + Likes + Dislikes + Comments + Aggregate.Followers, kernel.function = "gaussian", rho = 1.562652), data = csm)
#  

## ---- eval=LOCAL, echo = FALSE-------------------------------------------
#  
#  load(file = "D:/PostDoc_Canada/KSPM/csm.fit3.rda")
#  load(file = "D:/PostDoc_Canada/KSPM/summary.csm.fit3.rda")
#  
#  # summary(csm.fit3, kernel.test = "none", global.test = "TRUE")
#  
#  print(summary.csm.fit3)
#  

## ---- eval=LOCAL---------------------------------------------------------
#  # new data frame for Ker1
#  newdata.Ker1 <- data.frame(Genre = c(1, 3, 8), Gross = c(5.0e+07, 50000, 10000),Budget = c(1.8e+08, 5.2e+05, 1.3e+03), Screens = c(3600, 210, 5050), Sequel = c(2, 1, 1))
#  
#  # new data frame for Ker2
#  newdata.Ker2 <- data.frame(Sentiment = c(1, 2, 10), Views = c(293021, 7206, 5692061), Likes = c(3698, 2047, 5025), Dislikes = c(768, 49, 305), Comments = c(336, 70, 150), Aggregate.Followers = c(4530000, 350000, 960000))

## ---- eval=LOCAL---------------------------------------------------------
#  new.predictions <- predict(csm.fit3, newdata.kernel = list(Ker1 = newdata.Ker1, Ker2 = newdata.Ker2), interval = "prediction")
#  new.predictions

## ---- eval=LOCAL---------------------------------------------------------
#  head(csm.fit3$fitted.value)

## ---- eval=LOCAL---------------------------------------------------------
#  pred <- predict(csm.fit3, interval = "confidence")
#  head(pred)

## ---- fig.height=4, fig.show='hold', fig.width=4, eval=LOCAL-------------
#  plot(csm$Ratings, pred$fit, xlim = c(2, 10), ylim = c(2, 10), xlab = "Observed ratings", ylab = "Predicted ratings", cex.lab = 1.5)
#  abline(a = 0, b = 1, col = "red", lty = 2)

## ---- eval = FALSE-------------------------------------------------------
#  # NOT RUN
#  csm.fit4 <- kspm(response = "Ratings", linear = NULL, kernel = ~Kernel(~ Sentiment + Views + Likes + Dislikes + Comments + Aggregate.Followers, kernel.function = "gaussian"), data = csm, level = 0, control = kspmControl(parallel = TRUE))

## ---- eval = FALSE-------------------------------------------------------
#  # NOT RUN
#  stepKSPM(csm.fit4, kernel.lower = ~1, kernel.upper = ~ Sentiment + Views + Likes + Dislikes + Comments + Aggregate.Followers, direction = "both", k = 2, kernel.param = "change", data = csm)

## ---- eval=LOCAL---------------------------------------------------------
#  data("energy")
#  head(energy)

## ---- fig.height=4, fig.show='hold', fig.width=12, eval=LOCAL------------
#  par(mfrow = c(1,2), mar = c(5,5,2,2))
#  # energy among all the measurements
#  plot(energy$power, type = "l", xlab = "Time", ylab = "Power", cex.lab = 1.5, xaxt = "n")
#  axis(1, at = 1 + 24 * (0:21), labels = unique(energy$date))
#  # examples from three days
#  plot(c(NA,energy[1:26, "power"]), type = "b", xlab = "Time", ylab = "Power", cex.lab = 1.5, xaxt = "n", col = "darkgreen", lwd = 2, cex = 0.8, xlim = c(-1, 30))
#  lines(energy[24:50, "power"], type = "b", col = "blue", lwd = 2, cex = 0.8, pch = 0)
#  lines(energy[48:74, "power"], type = "b", col = "red", lwd = 2, cex = 1, pch = 17)
#  axis(1, at = c(1, 7, 13, 19, 25), labels = c("0h00", "6h00", "12h00", "18h00", "0h00"))
#  legend("topleft", col = c("darkgreen", "blue", "red"), legend = c("Sep 13, 2015", "Sep 14, 2015", "Sep 15, 2015"), lwd = 2, pch = c(1, 0, 17))
#  abline(v = 24.9, lty = 2)
#  text(25.5, 730, "next day", adj = 0)

## ---- eval=LOCAL---------------------------------------------------------
#  energy_train_ <- energy[1:408, ]
#  energy_test_ <- energy[409:504, ]

## ---- eval=LOCAL---------------------------------------------------------
#  energy.fit1 <- kspm(response = "power", linear = ~Temperature, kernel = ~Kernel(~hour.num + P + HR, kernel.function = "gaussian", rho = 0.7) , data = energy_train_)

## ---- eval=LOCAL---------------------------------------------------------
#  energy.fit1$kernel.info$Ker1$rho

## ---- eval=LOCAL---------------------------------------------------------
#  energy.fit2 <- kspm(response = "power", linear = ~Temperature, kernel = ~Kernel(~hour.num + P + HR, kernel.function = "gaussian", rho = 0.07) , data = energy_train_, level = 0)
#  energy.fit3 <- kspm(response = "power", linear = ~Temperature, kernel = ~Kernel(~hour.num + P + HR, kernel.function = "gaussian", rho = 7) , data = energy_train_, level = 0)

## ---- fig.height=15, fig.show='hold', fig.width=15, eval=LOCAL-----------
#  
#  ### parameters for figures panel
#  par(oma = c(1, 4, 6, 1))
#  par(mfrow = c(4,3), mar = c(5,5,1,1))
#  
#  ### kspm.fit1 (rho = 0.7)
#  # predictions with confidence intervals on train_
#  plot(energy_train_[1:72, "power"], type = "l", xlab = "Time", ylab = "Power", cex.lab = 1.5, xaxt = "n", lwd = 2, ylim = c(300, 750))
#  axis(1, at = 1 + 24 * (0:2), labels = unique(energy_train_$date)[1:3])
#  pred_train_ <- predict(energy.fit1, interval = "confidence")
#  lines(pred_train_$fit, col = "red")
#  lines(pred_train_$lwr, col = "blue", lty = 2)
#  lines(pred_train_$upr, col = "blue", lty = 2)
#  # predictions with prediction intervals on test_
#  plot(energy_test_[1:72, "power"], type = "l", xlab = "Time", ylab = "Power", cex.lab = 1.5, xaxt = "n", lwd = 2, ylim = c(300, 750))
#  axis(1, at = 1 + 24 * (0:2), labels = unique(energy_test_$date)[1:3])
#  pred_train_ <- predict(energy.fit1, newdata.linear =  energy_test_, newdata.kernel = list(Ker1 = energy_test_), interval = "prediction")
#  lines(pred_train_$fit, col = "red")
#  lines(pred_train_$lwr, col = "blue", lty = 2)
#  lines(pred_train_$upr, col = "blue", lty = 2)
#  # derivatives
#  plot(derivatives(energy.fit1), subset = "hour.num", xaxt = "n", ylab = "Derivatives", cex.lab = 1.5, ylim = c(-1000,1000))
#  axis(1, at = c(0, 6, 12, 18), labels = c("0h00", "6h00", "12h00", "18h00"))
#  
#  ### kspm.fit3 (rho = 0.07)
#  # predictions with confidence intervals on train_
#  plot(energy_train_[1:72, "power"], type = "l", xlab = "Time", ylab = "Power", cex.lab = 1.5, xaxt = "n", lwd = 2, ylim = c(300, 750))
#  axis(1, at = 1 + 24 * (0:2), labels = unique(energy_train_$date)[1:3])
#  pred_train_ <- predict(energy.fit2, interval = "confidence")
#  lines(pred_train_$fit, col = "red")
#  lines(pred_train_$lwr, col = "blue", lty = 2)
#  lines(pred_train_$upr, col = "blue", lty = 2)
#  # predictions with prediction intervals on test_
#  plot(energy_test_[1:72, "power"], type = "l", xlab = "Time", ylab = "Power", cex.lab = 1.5, xaxt = "n", lwd = 2, ylim = c(300, 750))
#  axis(1, at = 1 + 24 * (0:2), labels = unique(energy_test_$date)[1:3])
#  pred_train_ <- predict(energy.fit2, newdata.linear =  energy_test_, newdata.kernel = list(Ker1 = energy_test_), interval = "prediction")
#  lines(pred_train_$fit, col = "red")
#  lines(pred_train_$lwr, col = "blue", lty = 2)
#  lines(pred_train_$upr, col = "blue", lty = 2)
#  # derivatives
#  plot(derivatives(energy.fit2), subset = "hour.num", xaxt = "n", ylab = "Derivatives", cex.lab = 1.5, ylim = c(-1000,1000))
#  axis(1, at = c(0, 6, 12, 18), labels = c("0h00", "6h00", "12h00", "18h00"))
#  
#  ### kspm.fit2 (rho = 7)
#  # predictions with confidence intervals on train_
#  plot(energy_train_[1:72, "power"], type = "l", xlab = "Time", ylab = "Power", cex.lab = 1.5, xaxt = "n", lwd = 2, ylim = c(300, 750))
#  axis(1, at = 1 + 24 * (0:2), labels = unique(energy_train_$date)[1:3])
#  pred_train_ <- predict(energy.fit3, interval = "confidence")
#  lines(pred_train_$fit, col = "red")
#  lines(pred_train_$lwr, col = "blue", lty = 2)
#  lines(pred_train_$upr, col = "blue", lty = 2)
#  # predictions with prediction intervals on test_
#  plot(energy_test_[1:72, "power"], type = "l", xlab = "Time", ylab = "Power", cex.lab = 1.5, xaxt = "n", lwd = 2, ylim = c(300, 750))
#  axis(1, at = 1 + 24 * (0:2), labels = unique(energy_test_$date)[1:3])
#  pred_train_ <- predict(energy.fit3, newdata.linear =  energy_test_, newdata.kernel = list(Ker1 = energy_test_), interval = "prediction")
#  lines(pred_train_$fit, col = "red")
#  lines(pred_train_$lwr, col = "blue", lty = 2)
#  lines(pred_train_$upr, col = "blue", lty = 2)
#  # derivatives
#  plot(derivatives(energy.fit3), subset = "hour.num", xaxt = "n", ylab = "Derivatives", cex.lab = 1.5, ylim = c(-1000,1000))
#  axis(1, at = c(0, 6, 12, 18), labels = c("0h00", "6h00", "12h00", "18h00"))
#  
#  # Legends
#  plot.new()
#  legend("topleft", lty = c(1,1,2), col = c("black", "red", "blue"), legend = c("True data", "Predictions", "Confidence intervals"), cex = 2, bty = "n")
#  plot.new()
#  legend("topleft", lty = c(1,1,2), col = c("black", "red", "blue"), legend = c("True data", "Predictions", "Prediction intervals"), cex = 2,  bty = "n")
#  plot.new()
#  legend("topleft", pch = 1, col = c("black"), legend = c("Pointwise derivatives \n 1 point = 1 measure"), cex = 2, bty = "n")
#  
#  
#  ### legends on the left
#  par(fig = c(0,0.05,0,0.25), oma = c(0,0,0,0), mar = c(0,0,0,0),new = TRUE)
#  plot.new()
#  text(0.1, 0.8, "Legend", srt = 90, cex = 2, adj = 0.5)
#  par(fig = c(0,0.05,0.26,0.5), oma = c(0,0,0,0), mar = c(0,0,0,0),new = TRUE)
#  plot.new()
#  text(0.1, 0.5, expression(paste(rho, " = 7")), srt = 90, cex = 2, adj = 0.5)
#  par(fig = c(0,0.05,0.5,0.72), oma = c(0,0,0,0), mar = c(0,0,0,0),new = TRUE)
#  plot.new()
#  text(0.1, 0.5, expression(paste(rho, " = 0.07")), srt = 90, cex = 2, adj = 0.5)
#  par(fig = c(0,0.05,0.72,0.97), oma = c(0,0,0,0), mar = c(0,0,0,0),new = TRUE)
#  plot.new()
#  text(0.1, 0.5, expression(paste(rho, " = 0.7")), srt = 90, cex = 2, adj = 0.5)
#  
#  ### legends on the top
#  par(fig = c(0.05,0.36,0.92,1), oma = c(0,0,0,0), mar = c(0,0,0,0),new = TRUE)
#  plot.new()
#  text(0.5, 0.5, "Predictions on training data set",  cex = 2, adj = 0.5)
#  par(fig = c(0.37,0.68,0.92,1), oma = c(0,0,0,0), mar = c(0,0,0,0),new = TRUE)
#  plot.new()
#  text(0.5, 0.5, "Predictions on test data set",  cex = 2, adj = 0.5)
#  par(fig = c(0.69,1,0.92,1), oma = c(0,0,0,0), mar = c(0,0,0,0),new = TRUE)
#  plot.new()
#  text(0.5, 0.5, "Derivatives for hour.num",  cex = 2, adj = 0.5)

