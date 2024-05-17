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
fed = load.fed23()
LW = load.LW()
## ##################################
## Preprocessing
## ##################################
fed.fda = fda.preprocess(fed)
LW.fda = fda.preprocess(LW)
LW.fda$raw.data = LW.fda$raw.data[,c(1,10,seq(20,360, by=20))]
LW.fda$observationgrid = c(1,10, seq(20,360, by=20))
## ##################################
## 3D plots
## ##################################
pdf("figure3.pdf", width=36, height=16, pointsize = 30)
par(mfrow=c(1,2))
dffm.3Dplot(fed.fda, domainlab = "Maturity (months)", outputlab = "Yield (percent)")
dffm.3Dplot(LW.fda, domainlab = "Maturity (months)", outputlab = "Yield (percent)")
dev.off()



