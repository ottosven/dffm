## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Figure 5.
## ####################################################################
## ####################################################################
library(dffm)
## ##################################
## Data
## ##################################
fed = load.fed()
JKV = load.JKV()
## ##################################
## Preprocessing
## ##################################
JKV.fpca = fpca.preprocess(JKV, method="naturalsplines")
fed.fpca = fpca.preprocess(fed, method="naturalsplines")
## ##################################
## Criterion
## ##################################
JKV.criterion = dffm.criterion(JKV.fpca)
fed.criterion = dffm.criterion(fed.fpca)
## ##################################
## 3D plots
## ##################################
pdf("figure5.pdf", width=36, height=16, pointsize = 30)
par(mfrow=c(1,2))
dffm.MSEplot(JKV.criterion, kmax=8, pmax=8, zlim = range(na.omit(JKV.criterion$MSE.matrix[1:8,1:8])) + c(0,-11.5))
dffm.MSEplot(fed.criterion, kmax=8, pmax=8, zlim = range(na.omit(fed.criterion$MSE.matrix[1:8,1:8])) + c(0,-45))
dev.off()
