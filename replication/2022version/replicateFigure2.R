## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
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
JKV.fpca = fpca.preprocess(JKV, method="naturalsplines")
fed.fpca = fpca.preprocess(fed, method="naturalsplines")
## ##################################
## 3D plots
## ##################################
pdf("figure2.pdf", width=36, height=16, pointsize = 30)
par(mfrow=c(1,2))
dffm.3Dplot(JKV.fpca, domainlab = "Maturity (months)", outputlab = "Yield (percent)")
dffm.3Dplot(fed.fpca, domainlab = "Maturity (months)", outputlab = "Yield (percent)")
dev.off()
