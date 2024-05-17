## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Table 3.
## ####################################################################
## ####################################################################
library(dffm)
## ##################################
## Data
## ##################################
fed = load.fed23()
LW = load.LW()
fed.train = window(fed, end=time(fed)[120])
LW.train = window(LW, end=time(LW)[120])
## ##################################
## Preprocessing
## ##################################
fed.fda = fda.preprocess(fed, workinggrid = seq(1,360,by=0.5))
LW.fda = fda.preprocess(LW, workinggrid = seq(1,360,by=0.5))
fed.train.fda = fda.preprocess(fed.train, workinggrid = seq(1,360,by=0.5))
LW.train.fda = fda.preprocess(LW.train, workinggrid = seq(1,360,by=0.5))
## ##################################
## Information criteria
## ##################################
fed.criteria = fts.criterion(fed.fda, K.max=8, p.max=8)
LW.criteria = fts.criterion(LW.fda, K.max=8, p.max=8)
fed.train.criteria = fts.criterion(fed.train.fda, K.max=8, p.max=8)
LW.train.criteria = fts.criterion(LW.train.fda, K.max=8, p.max=8)
## ##################################
## Information criteria results
## ##################################
factorsandlags = rbind(
  c(fed.criteria$IC.min[-3,]),
  c(fed.train.criteria$IC.min[-3,]),
  c(LW.criteria$IC.min[-3,]),
  c(LW.train.criteria$IC.min[-3,])
)

rownames(factorsandlags) = c(
  "fed (full)",
  "fed (first 120)",
  "LW (full)",
  "LW (first 120)"
)
colnames(factorsandlags) = c(
  "K-bic", "K-hqc",
  "p-bic", "p-hqc"
)
t(factorsandlags)
write.table(t(factorsandlags),file=paste("./table3.csv", sep=''),append=T,col.names=F,row.names=F)
