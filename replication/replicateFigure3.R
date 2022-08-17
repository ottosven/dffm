## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
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
## Nelson-Siegel loading functions
## ##################################
figure.DNS = function(lambda = NULL, data.length = NULL, maturity.length = NULL){
  if(is.null(maturity.length)) maturity.length = 120
  if(is.null(data.length)) data.length = 360
  if(is.null(lambda)) lambda = 0.0609
  NSline <- matrix(nrow=data.length, ncol=3)
  for (i in 1:data.length){
    NSline[i, 1] <- 1
    NSline[i, 2] <- (1-exp(-lambda*i))/(lambda * i)
    NSline[i, 3] <- (1-exp(-lambda*i))/(lambda * i) - exp(- lambda * i)
  }
  plot(1:maturity.length, NSline[1:maturity.length,1], cex.lab=1.3, cex.axis=1.3, cex.main=1.6,
       type='l', ylim = c(0,1.5), lty=1, lwd = 4, xlab = "Maturity (months)", ylab = "", main= "Nelson-Siegel loading functions")
  lines(1:maturity.length, NSline[1:maturity.length,2], lty=2, lwd=4)
  lines(1:maturity.length, NSline[1:maturity.length,3], lty=4, lwd=3)
  legend(x = "top", legend = c('Level', 'Slope', 'Curvature'), lty = c(1,2,4),
         lwd = c(4,4,3), seg.len=2.5, cex = 1, pt.cex=2, horiz = TRUE, bty="n")
}
## ##################################
## FFM loading functions
## ##################################
figure.loadingfunctions = function(data, signpattern = c(1,1,1,1), main, ylim){
  FPCdata = fpca.preprocess(data, method="naturalsplines", workinggrid = seq(as.numeric(colnames(data))[1],as.numeric(colnames(data))[dim(data)[2]],0.5))
  FPCloadings = FPCdata$eigenfunctions.workgrid
  plot(FPCdata$workinggrid, signpattern[1]*FPCloadings[,1], ylim = ylim, type='l', lty=1, col=1, lwd=4,
       xlab = 'Maturity (months)', ylab = '', main=main,  cex.lab=1.3, cex.axis=1.3, cex.main=1.6)
  lines(FPCdata$workinggrid, signpattern[2]*FPCloadings[,2], lty=2, col=1, lwd = 4)
  lines(FPCdata$workinggrid, signpattern[3]*FPCloadings[,3], lty=4, col=1, lwd = 3)
  lines(FPCdata$workinggrid, signpattern[4]*FPCloadings[,4], lty=1, col="grey", lwd = 1)
  legend(x="top",
         legend = c(expression(hat(psi)[1]),expression(hat(psi)[2]), expression(hat(psi)[3]), expression(hat(psi)[4])),
         col = c(1,1,1,"grey"),
         lty = c(1,2,4,1),
         lwd = c(4,4,3,1),
         pt.cex = 2,
         cex = 1,
         text.col = "black",
         seg.len = 2.5,
         horiz = TRUE,
         bty = "n"
  )
}
## ##################################
## Loading functions DNS, JKV, FED
## ##################################
pdf("figure3.pdf", width=30, height=9, pointsize = 30)
par(mfrow=c(1,3))
figure.DNS()
figure.loadingfunctions(JKV, signpattern=c(1, 1, -1, -1), main="Estimated loading functions JKV data (K=4)", ylim = range(-0.22, 0.22))
figure.loadingfunctions(fed, signpattern=c(1, 1, -1, -1), main="Estimated loading functions FED data (K=4)", ylim = range(-0.15, 0.15))
dev.off()
