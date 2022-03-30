## ####################################################################
## ####################################################################
## Supplement for
## "Dynamic Factor Model for Functional Time Series:
## Identification, Estimation, and Prediction"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Figure 1.
## ####################################################################
## ####################################################################
library(dffm)
FED = load.fed()
crit.test = dffm.criterion(fpca.preprocess(diff(FED)))
crit.test$MSE.matrix = crit.test$MSE.matrix*3-0.1
crit.test$MSE.matrix[,7:8] = crit.test$MSE.matrix[,7:8] + 0.003
crit.test$MSE.matrix[8,7:8] = crit.test$MSE.matrix[8,7:8] + 0.003
crit.test$MSE.matrix[7,7:8] = crit.test$MSE.matrix[7,7:8] - 0.003
dffm.MSEplot(crit.test, kmax=8, pmax=8)
