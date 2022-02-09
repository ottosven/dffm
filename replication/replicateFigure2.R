## ####################################################################
## ####################################################################
## Supplement for
## "Dynamic Factor Model for Functional Time Series:
## Identification, Estimation, and Prediction"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Figure 2.
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
## 3D plots
## ##################################
dffm.3Dplot(JKV.fpca, domainlab = "Maturity (months)", outputlab = "Yield (percent)")
dffm.3Dplot(fed.fpca, domainlab = "Maturity (months)", outputlab = "Yield (percent)")
