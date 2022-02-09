## ####################################################################
## ####################################################################
## Supplement for
## "Dynamic Factor Model for Functional Time Series:
## Identification, Estimation, and Prediction"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Figure 3.
## ####################################################################
## ####################################################################
library(dffm)
## ##################################
## Data
## ##################################
fed = load.fed()
JKV = load.JKV()
## ##################################
## Nelson−Siegel loading functions
## ##################################
NSline <- matrix(nrow=360, ncol=3)
for (i in 1:360){
  NSline[i, 1] <- 1
  NSline[i, 2] <- (1-exp(-0.0609*i))/(0.0609 * i)
  NSline[i, 3] <- (1-exp(-0.0609*i))/(0.0609 * i) - exp(- 0.0609 * i)
}
plot(1:120, NSline[1:120,1], type='l', ylim = c(0,1.5), lty=1,
     lwd = 3, xlab = "Maturity (months)", ylab = "", main= "Nelson-Siegel loading functions")
lines(1:120, NSline[1:120,2], lty=2, lwd=3)
lines(1:120, NSline[1:120,3], lty=4, lwd=2)
legend(x = "top", legend = c('Level', 'Slope', 'Curvature'), lty = c(1,2,4),
       lwd = c(2,2,1), seg.len=2.5, cex = 1, pt.cex=2, horiz = TRUE, bty="n")
## ##################################
## Figure loading functions
## ##################################
figure.loadingfunctions = function(data, ...){
  FPCloadings = fpca.preprocess(data)$eigenfunctions.workgrid
  observationgrid = colnames(data)
  start = head(observationgrid,1)
  end = tail(observationgrid,1)
  plot(start:end, FPCloadings[,1], ylim = range(-1.5, 2.16), type='l', lty=1, col=1, lwd=3,
       xlab = 'Maturity (months)', ylab = '', ...)
  lines(start:end, FPCloadings[,2], lty=2, col=1, lwd = 3)
  lines(start:end, -FPCloadings[,3], lty=4, col=1, lwd = 2)
  lines(start:end, FPCloadings[,4], lty=1, col="grey", lwd = 1)
  legend(x="top",
         legend = c(expression(hat(psi)[1]),expression(hat(psi)[2]), expression(hat(psi)[3]), expression(hat(psi)[4])),
         col = c(1,1,1,"grey"),
         lty = c(1,2,4,1),
         lwd = c(3,3,2,1),
         pt.cex = 2,
         cex = 1,
         text.col = "black",
         seg.len = 2.5,
         horiz = TRUE,
         bty = "n"
  )
}
## ##################################
## Loading functions JKV
## ##################################
figure.loadingfunctions(JKV, main="Estimated loading functions JKV data (K=4)")
## ##################################
## Loading functions fed
## ##################################
figure.loadingfunctions(fed, main="Estimated loading functions FED data (K=4)")
