## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Figure 6.
## ####################################################################
## ####################################################################
library(dffm)
## ##################################
## Data
## ##################################
fed = load.fed23()
LW = load.LW()
## ##################################
## Preprocessing
## ##################################
LW.fda = fda.preprocess(LW)
fed.fda = fda.preprocess(fed)
## ##################################
## Criterion
## ##################################
LW.criterion = fts.criterion(LW.fda)
fed.criterion = fts.criterion(fed.fda)
## ##################################
## 3D plots
## ##################################
pdf("figure6.pdf", width=36, height=16, pointsize = 30)
par(mfrow=c(1,2))
dffm.MSEplot(fed.criterion, kmax=8, pmax=6, zlim = range(na.omit(fed.criterion$MSE.matrix[1:8,1:8])) + c(0,-57))
dffm.MSEplot(LW.criterion, kmax=8, pmax=6, zlim = range(na.omit(LW.criterion$MSE.matrix[1:8,1:8])) + c(0,-57))
dev.off()

