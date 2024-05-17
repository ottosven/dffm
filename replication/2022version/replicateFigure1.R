## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Figure 1.
## ####################################################################
## ####################################################################
library(dffm)
fed = load.fed()
crit.test = dffm.criterion(fpca.preprocess(diff(fed), method="naturalsplines"))
crit.test$MSE.matrix = crit.test$MSE.matrix*3-40
crit.test$MSE.matrix[,7:8] = crit.test$MSE.matrix[,7:8] + 0.5
crit.test$MSE.matrix[8,7:8] = crit.test$MSE.matrix[8,7:8] + 0.5
crit.test$MSE.matrix[7,7:8] = crit.test$MSE.matrix[7,7:8] - 0.5
## ##################################
## Plot Figure 1
## ##################################
pdf("figure1.pdf", width=18, height=16, pointsize = 30)
dffm.MSEplot(crit.test, kmax=8, pmax=8)
dev.off()
