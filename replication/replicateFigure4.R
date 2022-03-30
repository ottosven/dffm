## ####################################################################
## ####################################################################
## Supplement for
## "Dynamic Factor Model for Functional Time Series:
## Identification, Estimation, and Prediction"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Figure 4.
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
JKV.fpca = fpca.preprocess(JKV)
fed.fpca = fpca.preprocess(fed)
## ##################################
## Criterion
## ##################################
JKV.criterion = dffm.criterion(JKV.fpca)
fed.criterion = dffm.criterion(fed.fpca)
## ##################################
## 3D plots
## ##################################
dffm.MSEplot(JKV.criterion, kmax=8, pmax=8, zlim = range(na.omit(JKV.criterion$MSE.matrix[1:8,1:8])) + c(0,-0.1))
dffm.MSEplot(fed.criterion, kmax=8, pmax=8, zlim = range(na.omit(fed.criterion$MSE.matrix[1:8,1:8])) + c(0,-0.125))

