## ####################################################################
## ####################################################################
## Supplement for
## "Dynamic Factor Model for Functional Time Series:
## Identification, Estimation, and Prediction"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Table 3.
## ####################################################################
## ####################################################################
library(dffm)
## ##################################
## Data
## ##################################
fed = load.fed()
JKV = load.JKV()
JKV.train = window(JKV, end=time(JKV)[120])
fed.train = window(fed, end=time(fed)[120])
## ##################################
## Preprocessing
## ##################################
JKV.fpca = fpca.preprocess(JKV)
fed.fpca = fpca.preprocess(fed)
JKV.train.fpca = fpca.preprocess(JKV.train)
fed.train.fpca = fpca.preprocess(fed.train)
## ##################################
## Information criteria
## ##################################
JKV.criteria = dffm.criterion(JKV.fpca, K.max=8, p.max=8)
fed.criteria = dffm.criterion(fed.fpca, K.max=8, p.max=8)
JKV.train.criteria = dffm.criterion(JKV.train.fpca, K.max=8, p.max=8)
fed.train.criteria = dffm.criterion(fed.train.fpca, K.max=8, p.max=8)
## ##################################
## Information criteria results
## ##################################
factorsandlags = rbind(
  c(JKV.criteria$IC.min),
  c(fed.criteria$IC.min),
  c(JKV.train.criteria$IC.min),
  c(fed.train.criteria$IC.min)
)
rownames(factorsandlags) = c(
  "JKV (full data)",
  "fed (full data)",
  "JKV (first 120 months)",
  "fed (first 120 months)"
)
colnames(factorsandlags) = c(
  "K-bic", "K-hqc", "K-fFPE",
  "p-bic", "p-hqc", "p-fFPE"
)
factorsandlags
