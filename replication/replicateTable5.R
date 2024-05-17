## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Table 5.
## THIS TAKES APPROXIMATELY 11 HOURS ON MY MACHINE
## The results are available in the folder "table5-results"
## ####################################################################
## ####################################################################
rm(list=ls())
start = Sys.time()
library(fdadata)

## ##################################
## Load Data
## ##################################
load("./not-for-public-access/eikonyields/G7data.RData")
FED = load.fed23()
LW = load.LW()
## #########################################################
## in-sample
## #########################################################
get.RMSFE.insamp = function(fdaobj, fit.VAR){
  p = dim(fit.VAR$VARmatrix)[2]/dim(fit.VAR$VARmatrix)[1]
  if(dim(fdaobj$raw.data)[2] == dim(fit.VAR$curve.predict)[2]){
    curve.pred = fit.VAR$curve.predict
  } else {
    curve.pred = fit.VAR$curve.predict[,match(fdaobj$observationgrid, fdaobj$workinggrid)]
  }
  sqrt(mean((curve.pred - fdaobj$raw.data[-(1:p),])^2, na.rm=TRUE))
}

get.rmsetable.insamp = function(data){
  workinggrid = seq(from=as.numeric(colnames(data)[1]),to=as.numeric(rev(colnames(data))[1]), by=0.5)
  fdaobj = fda.preprocess(data, workinggrid = workinggrid)
  cumAC = fts.cumAC(fdaobj)
  IC = fts.criterion(fdaobj, K.max = 12, p.max = 12)$IC.min
  rmse.BICVAR = get.RMSFE.insamp(fdaobj, fts.VARforecast(
    cumAC, K = IC[1,1], p = IC[1,2], AR = FALSE, start = NULL, end = NULL, h=0)
  )
  rmse.BICAR = get.RMSFE.insamp(fdaobj, fts.VARforecast(
    cumAC, K = IC[1,1], p = IC[1,2], AR = TRUE, start = NULL, end = NULL, h=0)
  )
  rmse.HQCVAR = get.RMSFE.insamp(fdaobj, fts.VARforecast(
    cumAC, K = IC[2,1], p = IC[2,2], AR = FALSE, start = NULL, end = NULL, h=0)
  )
  rmse.HQCAR = get.RMSFE.insamp(fdaobj, fts.VARforecast(
    cumAC, K = IC[2,1], p = IC[2,2], AR = TRUE, start = NULL, end = NULL, h=0)
  )
  rmse.AC3ARp1 = get.RMSFE.insamp(fdaobj, fts.VARforecast(
    cumAC, K = 3, p = 1, AR = TRUE, start = NULL, end = NULL, h=0)
  )
  rmse.AC4ARp1 = get.RMSFE.insamp(fdaobj, fts.VARforecast(
    cumAC, K = 4, p = 1, AR = TRUE, start = NULL, end = NULL, h=0)
  )
  rmse.AC6ARp1 = get.RMSFE.insamp(fdaobj, fts.VARforecast(
    cumAC, K = 6, p = 1, AR = TRUE, start = NULL, end = NULL, h=0)
  )
  rmse.AC8ARp1 = get.RMSFE.insamp(fdaobj, fts.VARforecast(
    cumAC, K = 8, p = 1, AR = TRUE, start = NULL, end = NULL, h=0)
  )
  rmse.DNSVARp1 = get.RMSFE.insamp(fdaobj, DNS.forecast(
    fdaobj, p = 1, AR = FALSE, obsdata = FALSE, start = NULL, end = NULL, h=0)
  )
  rmse.DNSARp1 = get.RMSFE.insamp(fdaobj, DNS.forecast(
    fdaobj, p = 1, AR = TRUE, obsdata = FALSE, start = NULL, end = NULL, h=0)
  )
  rmse.RW = sqrt(mean((diff(data,1))^2, na.rm=TRUE))
  ### table
  rmse.all = round(c(
    rmse.AC3ARp1,
    rmse.AC4ARp1,
    rmse.AC6ARp1,
    rmse.AC8ARp1,
    rmse.BICAR,
    rmse.HQCAR,
    rmse.BICVAR,
    rmse.HQCVAR,
    rmse.DNSVARp1,
    rmse.DNSARp1,
    rmse.RW
  ),3)
  names(rmse.all) = c(
    "AR3", "AR4", "AR6", "AR8", "BICAR", "HQCAR",
    "BICVAR", "HQCVAR",
    "DNSVAR", "DNSAR",
    "RW"
  )
  return(rmse.all)
}

## #########################################################
## h-step predictions
## #########################################################

get.MSFE.oos = function(i, fdaobj, fdapred, K=3, AR=FALSE, h=1){
  pred.VAR = fts.VARforecast(fdapred, K = K, p = 1, AR = AR, end = i, h=h)
  curve.pred = pred.VAR[match(fdaobj$observationgrid, fdaobj$workinggrid)]
  return(mean((curve.pred - fdaobj$raw.data[i+h,])^2, na.rm=TRUE))
}

get.MSFE.oos.DNS = function(i, fdaobj, AR=FALSE, h=1){
  pred.VAR = DNS.forecast(fdaobj, p=1, AR=AR, end=i, h=h)
  curve.pred = pred.VAR[match(fdaobj$observationgrid, fdaobj$workinggrid)]
  return(mean((curve.pred - fdaobj$raw.data[i+h,])^2, na.rm=TRUE))
}

get.MSFE.oos.IC = function(i, fdaobj, fdapred, h=1){
  re.fdaobj = fda.preprocess(
    window(fdaobj$raw.data, end = time(fdaobj$raw.data)[i]),
    workinggrid = fdaobj$workinggrid)
  IC = fts.criterion(re.fdaobj, K.max = 12, p.max = 12)$IC.min
  pred.BICVAR = fts.VARforecast(fdapred, K = IC[1,1], p = IC[1,2], AR = FALSE, end = i, h=h)
  curve.pred.BICVAR = pred.BICVAR[match(fdaobj$observationgrid, fdaobj$workinggrid)]
  pred.BICAR = fts.VARforecast(fdapred, K = IC[1,1], p = IC[1,2], AR = TRUE, end = i, h=h)
  curve.pred.BICAR = pred.BICAR[match(fdaobj$observationgrid, fdaobj$workinggrid)]
  pred.HQCVAR = fts.VARforecast(fdapred, K = IC[2,1], p = IC[2,2], AR = FALSE, end = i, h=h)
  curve.pred.HQCVAR = pred.HQCVAR[match(fdaobj$observationgrid, fdaobj$workinggrid)]
  pred.HQCAR = fts.VARforecast(fdapred, K = IC[2,1], p = IC[2,2], AR = TRUE, end = i, h=h)
  curve.pred.HQCAR = pred.HQCAR[match(fdaobj$observationgrid, fdaobj$workinggrid)]
  mse.out = c(
    mean((curve.pred.BICAR - fdaobj$raw.data[i+h,])^2, na.rm=TRUE),
    mean((curve.pred.HQCAR - fdaobj$raw.data[i+h,])^2, na.rm=TRUE),
    mean((curve.pred.BICVAR - fdaobj$raw.data[i+h,])^2, na.rm=TRUE),
    mean((curve.pred.HQCVAR - fdaobj$raw.data[i+h,])^2, na.rm=TRUE)
  )
  return(mse.out)
}

get.rmsetable.oos = function(h, data, training = 120){
  workinggrid = seq(from=as.numeric(colnames(data)[1]),to=as.numeric(rev(colnames(data))[1]), by=0.5)
  fdaobj = fda.preprocess(data, workinggrid = workinggrid)
  cumAC = fts.cumAC(fdaobj)
  lastobs = dim(fdaobj$raw.data)[1]-h
  period = training:(dim(fdaobj$raw.data)[1]-h)
  rmses.autoselect = sqrt(rowMeans(sapply(period, get.MSFE.oos.IC, fdaobj=fdaobj, fdapred = cumAC, h=h)))
  rmse.AC3ARp1 = sqrt(mean(
    sapply(period, get.MSFE.oos, fdaobj=fdaobj, fdapred = cumAC, K=3, AR=TRUE, h=h)))
  rmse.AC4ARp1 = sqrt(mean(
    sapply(period, get.MSFE.oos, fdaobj=fdaobj, fdapred = cumAC, K=4, AR=TRUE, h=h)))
  rmse.AC6ARp1 = sqrt(mean(
    sapply(period, get.MSFE.oos, fdaobj=fdaobj, fdapred = cumAC, K=6, AR=TRUE, h=h)))
  rmse.AC8ARp1 = sqrt(mean(
    sapply(period, get.MSFE.oos, fdaobj=fdaobj, fdapred = cumAC, K=8, AR=TRUE, h=h)))
  rmse.DNSVARp1 = sqrt(mean(
    sapply(period, get.MSFE.oos.DNS, fdaobj=fdaobj, AR=FALSE, h=h)))
  rmse.DNSARp1 = sqrt(mean(
    sapply(period, get.MSFE.oos.DNS, fdaobj=fdaobj, AR=TRUE, h=h)))
  rmse.RW = sqrt(mean((diff(data,h))^2, na.rm=TRUE))
  ### table
  rmse.all = round(c(
    rmse.AC3ARp1,
    rmse.AC4ARp1,
    rmse.AC6ARp1,
    rmse.AC8ARp1,
    rmses.autoselect,
    rmse.DNSVARp1,
    rmse.DNSARp1,
    rmse.RW
  ),3)
  names(rmse.all) = c(
    "AR3", "AR4", "AR6", "AR8", "BICAR", "HQCAR",
    "BICVAR", "HQCVAR",
    "DNSVAR", "DNSAR",
    "RW"
  )
  rmse.all
  return(rmse.all)
}


## #########################################################
## Save results in "table5-results"
## #########################################################

data = G7data$FR
insamp = get.rmsetable.insamp(data)
oos = sapply(c(1,3,6,12), get.rmsetable.oos, data=data)
results = t(cbind(insamp, oos))
results
rownames(results) = c("insample", "1step", "3step", "6step", "12step")
results
write.csv(results, file = "./table5-results/results-FR.csv")
Sys.time() - start

data = G7data$DE
insamp = get.rmsetable.insamp(data)
oos = sapply(c(1,3,6,12), get.rmsetable.oos, data=data)
results = t(cbind(insamp, oos))
rownames(results) = c("insample", "1step", "3step", "6step", "12step")
results
write.csv(results, file = "./table5-results/results-DE.csv")
Sys.time() - start

data = G7data$IT
insamp = get.rmsetable.insamp(data)
oos = sapply(c(1,3,6,12), get.rmsetable.oos, data=data)
results = t(cbind(insamp, oos))
rownames(results) = c("insample", "1step", "3step", "6step", "12step")
results
write.csv(results, file = "./table5-results/results-IT.csv")
Sys.time() - start

data = G7data$JP
insamp = get.rmsetable.insamp(data)
oos = sapply(c(1,3,6,12), get.rmsetable.oos, data=data)
results = t(cbind(insamp, oos))
rownames(results) = c("insample", "1step", "3step", "6step", "12step")
results
write.csv(results, file = "./table5-results/results-JP.csv")
Sys.time() - start

data = G7data$GB
insamp = get.rmsetable.insamp(data)
oos = sapply(c(1,3,6,12), get.rmsetable.oos, data=data)
results = t(cbind(insamp, oos))
rownames(results) = c("insample", "1step", "3step", "6step", "12step")
results
write.csv(results, file = "./table5-results/results-GB.csv")
Sys.time() - start

data = G7data$US
insamp = get.rmsetable.insamp(data)
oos = sapply(c(1,3,6,12), get.rmsetable.oos, data=data)
results = t(cbind(insamp, oos))
rownames(results) = c("insample", "1step", "3step", "6step", "12step")
results
write.csv(results, file = "./table5-results/results-US.csv")
Sys.time() - start

data = G7data$CA
insamp = get.rmsetable.insamp(data)
oos = sapply(c(1,3,6,12), get.rmsetable.oos, data=data)
results = t(cbind(insamp, oos))
rownames(results) = c("insample", "1step", "3step", "6step", "12step")
results
write.csv(results, file = "./table5-results/results-CA.csv")
Sys.time() - start

data = FED
insamp = get.rmsetable.insamp(data)
oos = sapply(c(1,3,6,12), get.rmsetable.oos, data=data)
results = t(cbind(insamp, oos))
rownames(results) = c("insample", "1step", "3step", "6step", "12step")
results
write.csv(results, file = "./table5-results/results-FED.csv")
Sys.time() - start

data = LW
insamp = get.rmsetable.insamp(data)
oos = sapply(c(1,3,6,12), get.rmsetable.oos, data=data)
results = t(cbind(insamp, oos))
rownames(results) = c("insample", "1step", "3step", "6step", "12step")
results
write.csv(results, file = "./table5-results/results-LW.csv")
Sys.time() - start
