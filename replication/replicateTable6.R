## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Table 6.
## ####################################################################
## ####################################################################
## NOTE: This code runs about 3-5 hours
## ####################################################################
## The G7 data used in this paper are the property of Thomson Reuters Eikon
## and must be purchased from Thomson Reuters Eikon to gain access.
## ####################################################################
## To reproduce the G7 dataset, one has to follow these steps:
## 1) You have access to a Thomson Reuters Eikon account
## 2) You have Thomson Reuters Eikon and the Eikon Excel Add-In installed
## 3) Use the Excel Add-In to retrieve the data. 
## Each time series in Thomson Reuters Eikon has a specific code. 
## For instance, CA1YT=RR are the yields of a zero coupon bond for Canada with 1 year time to maturity
## Use the following codes to retrieve the time series:
## Canada: CA1MT=RR, CA2MT=RR, CA3MT=RR, CA6MT=RR, CA1YT=RR, CA2YT=RR, CA3YT=RR, CA4YT=RR, CA5YT=RR, CA7YT=RR, CA10YT=RR, CA20YT=RR, CA30YT=RR
## France: FR2YT=RR, FR3YT=RR, FR4YT=RR, FR5YT=RR, FR6YT=RR, FR7YT=RR, FR8YT=RR, FR9YT=RR, FR10YT=RR, FR15YT=RR, FR20YT=RR, FR25YT=RR, FR30YT=RR
## Germany: DE2YT=RR, DE3YT=RR, DE4YT=RR, DE5YT=RR, DE6YT=RR, DE7YT=RR, DE8YT=RR, DE9YT=RR, DE10YT=RR, DE15YT=RR, DE20YT=RR, DE25YT=RR, DE30YT=RR
## Italy: IT3YT=RR, IT4YT=RR, IT5YT=RR, IT6YT=RR, IT7YT=RR, IT8YT=RR, IT9YT=RR, IT10YT=RR, IT15YT=RR, IT20YT=RR, IT25YT=RR, IT30YT=RR
## Japan: JP3MT=RR, JP6MT=RR, JP9MT=RR, JP1YT=RR, JP2YT=RR, JP3YT=RR, JP4YT=RR, JP5YT=RR, JP6YT=RR, JP7YT=RR, JP8YT=RR, JP9YT=RR, JP10YT=RR, JP15YT=RR, JP20YT=RR
## United Kingdom: GB3MT=RR, GB6MT=RR, GB1YT=RR, GB2YT=RR, GB3YT=RR, GB4YT=RR, GB5YT=RR, GB6YT=RR, GB7YT=RR, GB8YT=RR, GB9YT=RR, GB10YT=RR, GB12YT=RR, GB15YT=RR, GB20YT=RR, GB25YT=RR, GB30YT=RR
## United States: US3MT=RR, US6MT=RR, US1YT=RR, US2YT=RR, US3YT=RR, US5YT=RR, US7YT=RR, US10YT=RR, US20YT=RR, US30YT=RR
## For each time series use the following parameters:
## Start:01.01.1995, End:30.06.2022, Interval:1MO
## 4) Save the data for each country into a mts-object with monthly frequency. The column names are the times to maturity in months. 
## Note that some values are missing in the original dataset. Missing vallues in the mts-object are set to "NA".
## 5) Save the list of the seven mts-objects as "G7data.Rdata". The names of the list elements are c("CA", "FR", "DE", "IT", "JP", "GB", "US")
## ####################################################################
## ####################################################################
library(dffm)
## ##################################
## Data
## ##################################
load("./not-for-public-access/G7data.RData")
## ##################################
## ##################################
## In sample forecast functions
## ##################################
dffm.insampleforecasts = function(dffmobj){
  factordynamics = match.arg(dffmobj$factordynamics, c("VAR","AR"))
  if(factordynamics == "AR"){
    VARfitted = matrix(nrow = dim(dffmobj$factors)[1]-dffmobj$p, ncol=dffmobj$K)
    for(i in 1:dffmobj$K) VARfitted[,i] = lm(embed(dffmobj$factors[,i],dffmobj$p+1)[,1] ~ -1 + embed(dffmobj$factors[,i],dffmobj$p+1)[,-1])$fitted.values
  } else {
    VARfitted = lm(embed(dffmobj$factors, dffmobj$p+1)[,1:dffmobj$K] ~ -1 + embed(dffmobj$factors, dffmobj$p+1)[,-(1:dffmobj$K)])$fitted.values
  }
  fittedcurve.obsgrid = VARfitted %*% t(dffmobj$loadingfunctions.obsgrid) + matrix(rep(dffmobj$meanfunction.obsgrid, dim(VARfitted)[1]), nrow = dim(VARfitted)[1], byrow=TRUE)
  return(fittedcurve.obsgrid)
}
get.workinggrid = function(data, gridsize){
  obsgrid = as.numeric(colnames(data))
  seq(obsgrid[1], obsgrid[length(obsgrid)], gridsize)
}
get.preprocess = function(data){
  fpca.preprocess(data, method = "naturalsplines", workinggrid = get.workinggrid(data, 0.5))
}
insample.prederror = function(data){
  prederrors = list()
  for(j in 1:12) prederrors[[j]] = matrix(nrow = dim(data)[1]-1, ncol = dim(data)[2])
  data.fpca = get.preprocess(data)
  IC = dffm.criterion(data.fpca)$IC.min
  prederrors[[1]] = dffm.insampleforecasts(dffm(data.fpca, K=IC[1,1], p=IC[1,2])) - data[-(1:IC[1,2]),]
  prederrors[[2]] = dffm.insampleforecasts(dffm(data.fpca, K=IC[2,1], p=IC[2,2])) - data[-(1:IC[2,2]),]
  prederrors[[3]] = dffm.insampleforecasts(dffm(data.fpca, K=IC[1,1], p=IC[1,2], AR=TRUE)) - data[-(1:IC[1,2]),]
  prederrors[[4]] = dffm.insampleforecasts(dffm(data.fpca, K=IC[2,1], p=IC[2,2], AR=TRUE)) - data[-(1:IC[2,2]),]
  prederrors[[5]] = dffm.insampleforecasts(dffm(data.fpca, K=3, p=1)) - data[-1,]
  prederrors[[6]] = dffm.insampleforecasts(dffm(data.fpca, K=4, p=1)) - data[-1,]
  prederrors[[7]] = dffm.insampleforecasts(dffm(data.fpca, K=6, p=1)) - data[-1,]
  prederrors[[8]] = dffm.insampleforecasts(dffm(data.fpca, K=3, p=1, AR=TRUE)) - data[-1,]
  prederrors[[9]] = dffm.insampleforecasts(dffm(data.fpca, K=4, p=1, AR=TRUE)) - data[-1,]
  prederrors[[10]] = dffm.insampleforecasts(dffm(data.fpca, K=6, p=1, AR=TRUE)) - data[-1,]
  prederrors[[11]] = DNS.insample(data, p=1, AR=FALSE) - data[-1,]
  prederrors[[12]] = DNS.insample(data, p=1, AR=TRUE) - data[-1,]
  rmse = numeric(12)
  for(j in 1:12) rmse[j] = sqrt(mean(prederrors[[j]]^2, na.rm=TRUE))
  return(rmse)
}
# Function for Nelson-Siegel loadings
DNS.insample = function(data, p=1, AR=FALSE){
  maturities = as.numeric(colnames(data))
  NS <- matrix(nrow=length(maturities), ncol=3)
  lambda = 0.0609
  for (i in 1:length(maturities)){
    NS[i, 1] <- 1
    NS[i, 2] <- (1-exp(-lambda*maturities[i]))/(lambda * maturities[i])
    NS[i, 3] <- (1-exp(-lambda*maturities[i]))/(lambda * maturities[i]) - exp(- lambda * maturities[i])}
  beta <- matrix(nrow=dim(data)[1], ncol = 3)
  for (t in 1:dim(data)[1]){
    beta[t,] <- lm(data[t,] ~ NS - 1)$coefficients
  }
  if(AR==TRUE){
    ARfit1 = lm(embed(beta[,1],p+1)[,1] ~ embed(beta[,1],p+1)[,-1])$fitted.values
    ARfit2 = lm(embed(beta[,2],p+1)[,1] ~ embed(beta[,2],p+1)[,-1])$fitted.values
    ARfit3 = lm(embed(beta[,3],p+1)[,1] ~ embed(beta[,3],p+1)[,-1])$fitted.values
    VARfitted = cbind(ARfit1, ARfit2, ARfit3)
  } else {
    VARfitted = lm(embed(beta, p+1)[,1:3] ~ embed(beta, p+1)[,-(1:3)])$fitted.values
  }
  fittedcurve = VARfitted %*% t(NS)
  return(fittedcurve)
}
## ##################################
## Out of sample forecast functions
## ##################################
outofsample.prederror = function(data, h=1, traininglength = 120){
  trainingendperiods = (traininglength):(dim(data)[1]-h)
  evalperiods = trainingendperiods + h
  prederrors = list()
  for(j in 1:12) prederrors[[j]] = matrix(nrow = length(trainingendperiods), ncol = dim(data)[2])
  for(i in 1:length(evalperiods)){
    traindata = window(data, end=time(data)[trainingendperiods[i]])
    data.fpca = get.preprocess(traindata)
    IC = dffm.criterion(data.fpca)$IC.min
    prederrors[[1]][i,] = tail(dffm.forecast(dffm(data.fpca, K=IC[1,1], p=IC[1,2]), h=h)$predcurves.obsgrid,1) - data[evalperiods[i],]
    prederrors[[2]][i,] = tail(dffm.forecast(dffm(data.fpca, K=IC[2,1], p=IC[2,2]), h=h)$predcurves.obsgrid,1) - data[evalperiods[i],]
    prederrors[[3]][i,] = tail(dffm.forecast(dffm(data.fpca, K=IC[1,1], p=IC[1,2], AR=TRUE), h=h)$predcurves.obsgrid,1) - data[evalperiods[i],]
    prederrors[[4]][i,] = tail(dffm.forecast(dffm(data.fpca, K=IC[2,1], p=IC[2,2], AR=TRUE), h=h)$predcurves.obsgrid,1) - data[evalperiods[i],]
    prederrors[[5]][i,] = tail(dffm.forecast(dffm(data.fpca, K=3, p=1), h=h)$predcurves.obsgrid,1) - data[evalperiods[i],]
    prederrors[[6]][i,] = tail(dffm.forecast(dffm(data.fpca, K=4, p=1), h=h)$predcurves.obsgrid,1) - data[evalperiods[i],]
    prederrors[[7]][i,] = tail(dffm.forecast(dffm(data.fpca, K=6, p=1), h=h)$predcurves.obsgrid,1) - data[evalperiods[i],]
    prederrors[[8]][i,] = tail(dffm.forecast(dffm(data.fpca, K=3, p=1, AR=TRUE), h=h)$predcurves.obsgrid,1) - data[evalperiods[i],]
    prederrors[[9]][i,] = tail(dffm.forecast(dffm(data.fpca, K=4, p=1, AR=TRUE), h=h)$predcurves.obsgrid,1) - data[evalperiods[i],]
    prederrors[[10]][i,] = tail(dffm.forecast(dffm(data.fpca, K=6, p=1, AR=TRUE), h=h)$predcurves.obsgrid,1) - data[evalperiods[i],]
    prederrors[[11]][i,] = DNS.forecast(traindata, p=1, h=h, AR=FALSE) - data[evalperiods[i],]
    prederrors[[12]][i,] = DNS.forecast(traindata, p=1, h=h, AR=TRUE) - data[evalperiods[i],]
  }
  rmse = numeric(12)
  for(j in 1:12) rmse[j] = sqrt(mean(prederrors[[j]]^2, na.rm=TRUE))
  return(rmse)
}
# Function for Nelson-Siegel loadings
DNS.forecast = function(data, p=1, h=1, AR=FALSE){
  maturities = as.numeric(colnames(data))
  NS <- matrix(nrow=length(maturities), ncol=3)
  lambda = 0.0609
  for (i in 1:length(maturities)){
    NS[i, 1] <- 1
    NS[i, 2] <- (1-exp(-lambda*maturities[i]))/(lambda * maturities[i])
    NS[i, 3] <- (1-exp(-lambda*maturities[i]))/(lambda * maturities[i]) - exp(- lambda * maturities[i])}
  beta <- matrix(nrow=dim(data)[1], ncol = 3)
  for (t in 1:dim(data)[1]){
    beta[t,] <- lm(data[t,] ~ NS - 1)$coefficients
  }
  if(AR==TRUE){
    VARcoefficients = matrix(0, nrow = 3*p+1, ncol=3)
    for(i in 1:3){
      ARfit = lm(embed(beta[,i],p+1)[,1] ~ embed(beta[,i],p+1)[,-1])
      VARcoefficients[1,i] = ARfit$coefficients[1]
      for(j in 1:p){
        VARcoefficients[i+3*(j-1)+1,i] = ARfit$coefficients[j+1]
      }
    }
  } else {
    VARcoefficients = lm(embed(beta, p+1)[,1:3] ~ embed(beta, p+1)[,-(1:3)])$coefficients
  }
  betahat = beta
  for(i in 1:h){
    betahat = rbind(betahat, c(1,tail(embed(betahat,p),1)) %*% VARcoefficients)
  }
  forecast = tail(betahat,1) %*% t(NS)
  return(forecast)
}
## ##################################
## Results for DE
## ##################################
start=Sys.time()
currentdata = G7data$DE
insample = insample.prederror(currentdata)
out1ahead = outofsample.prederror(currentdata, h=1, traininglength=120)
out3head = outofsample.prederror(currentdata, h=3, traininglength=120)
out6ahead = outofsample.prederror(currentdata, h=6, traininglength=120)
RMSFE = round(rbind(insample, out1ahead, out3head, out6ahead),3)
nnames=c(
  "VAR-BIC", "VAR-HQC", "AR-BIC", "AR-HQC",
  paste("VAR_",c(3,4,6),"(1)",sep=""),
  paste("AR_",c(3,4,6),"(1)",sep=""),
  "DNS-VAR", "DNS-AR"
)
colnames(RMSFE) = nnames
rownames(RMSFE) = c("in-sample 1-step", "1-ahead", "3-ahead", "6-ahead")
RMSFE.JKV
Sys.time()-start
write.table(RMSFE,file=paste("./table6-DE.csv", sep=''),append=T,col.names=F,row.names=F)
## ##################################
## Results for CA
## ##################################
start=Sys.time()
currentdata = G7data$CA
insample = insample.prederror(currentdata)
out1ahead = outofsample.prederror(currentdata, h=1, traininglength=120)
out3head = outofsample.prederror(currentdata, h=3, traininglength=120)
out6ahead = outofsample.prederror(currentdata, h=6, traininglength=120)
RMSFE = round(rbind(insample, out1ahead, out3head, out6ahead),3)
nnames=c(
  "VAR-BIC", "VAR-HQC", "AR-BIC", "AR-HQC",
  paste("VAR_",c(3,4,6),"(1)",sep=""),
  paste("AR_",c(3,4,6),"(1)",sep=""),
  "DNS-VAR", "DNS-AR"
)
colnames(RMSFE) = nnames
rownames(RMSFE) = c("in-sample 1-step", "1-ahead", "3-ahead", "6-ahead")
RMSFE
Sys.time()-start
write.table(RMSFE,file=paste("./table6-CA.csv", sep=''),append=T,col.names=F,row.names=F)
## ##################################
## Results for FR
## ##################################
start=Sys.time()
currentdata = G7data$FR
insample = insample.prederror(currentdata)
out1ahead = outofsample.prederror(currentdata, h=1, traininglength=120)
out3head = outofsample.prederror(currentdata, h=3, traininglength=120)
out6ahead = outofsample.prederror(currentdata, h=6, traininglength=120)
RMSFE = round(rbind(insample, out1ahead, out3head, out6ahead),3)
nnames=c(
  "VAR-BIC", "VAR-HQC", "AR-BIC", "AR-HQC",
  paste("VAR_",c(3,4,6),"(1)",sep=""),
  paste("AR_",c(3,4,6),"(1)",sep=""),
  "DNS-VAR", "DNS-AR"
)
colnames(RMSFE) = nnames
rownames(RMSFE) = c("in-sample 1-step", "1-ahead", "3-ahead", "6-ahead")
RMSFE
Sys.time()-start
write.table(RMSFE,file=paste("./table6-FR.csv", sep=''),append=T,col.names=F,row.names=F)
## ##################################
## Results for IT
## ##################################
start=Sys.time()
currentdata = G7data$IT
insample = insample.prederror(currentdata)
out1ahead = outofsample.prederror(currentdata, h=1, traininglength=120)
out3head = outofsample.prederror(currentdata, h=3, traininglength=120)
out6ahead = outofsample.prederror(currentdata, h=6, traininglength=120)
RMSFE = round(rbind(insample, out1ahead, out3head, out6ahead),3)
nnames=c(
  "VAR-BIC", "VAR-HQC", "AR-BIC", "AR-HQC",
  paste("VAR_",c(3,4,6),"(1)",sep=""),
  paste("AR_",c(3,4,6),"(1)",sep=""),
  "DNS-VAR", "DNS-AR"
)
colnames(RMSFE) = nnames
rownames(RMSFE) = c("in-sample 1-step", "1-ahead", "3-ahead", "6-ahead")
RMSFE
Sys.time()-start
write.table(RMSFE,file=paste("./table6-IT.csv", sep=''),append=T,col.names=F,row.names=F)
## ##################################
## Results for JP
## ##################################
start=Sys.time()
currentdata = G7data$JP
insample = insample.prederror(currentdata)
out1ahead = outofsample.prederror(currentdata, h=1, traininglength=120)
out3head = outofsample.prederror(currentdata, h=3, traininglength=120)
out6ahead = outofsample.prederror(currentdata, h=6, traininglength=120)
RMSFE = round(rbind(insample, out1ahead, out3head, out6ahead),3)
nnames=c(
  "VAR-BIC", "VAR-HQC", "AR-BIC", "AR-HQC",
  paste("VAR_",c(3,4,6),"(1)",sep=""),
  paste("AR_",c(3,4,6),"(1)",sep=""),
  "DNS-VAR", "DNS-AR"
)
colnames(RMSFE) = nnames
rownames(RMSFE) = c("in-sample 1-step", "1-ahead", "3-ahead", "6-ahead")
RMSFE
Sys.time()-start
write.table(RMSFE,file=paste("./table6-JP.csv", sep=''),append=T,col.names=F,row.names=F)
## ##################################
## Results for GB
## ##################################
start=Sys.time()
currentdata = G7data$GB
insample = insample.prederror(currentdata)
out1ahead = outofsample.prederror(currentdata, h=1, traininglength=120)
out3head = outofsample.prederror(currentdata, h=3, traininglength=120)
out6ahead = outofsample.prederror(currentdata, h=6, traininglength=120)
RMSFE = round(rbind(insample, out1ahead, out3head, out6ahead),3)
nnames=c(
  "VAR-BIC", "VAR-HQC", "AR-BIC", "AR-HQC",
  paste("VAR_",c(3,4,6),"(1)",sep=""),
  paste("AR_",c(3,4,6),"(1)",sep=""),
  "DNS-VAR", "DNS-AR"
)
colnames(RMSFE) = nnames
rownames(RMSFE) = c("in-sample 1-step", "1-ahead", "3-ahead", "6-ahead")
RMSFE
Sys.time()-start
write.table(RMSFE,file=paste("./table6-GB.csv", sep=''),append=T,col.names=F,row.names=F)
## ##################################
## Results for US
## ##################################
start=Sys.time()
currentdata = G7data$US
insample = insample.prederror(currentdata)
out1ahead = outofsample.prederror(currentdata, h=1, traininglength=120)
out3head = outofsample.prederror(currentdata, h=3, traininglength=120)
out6ahead = outofsample.prederror(currentdata, h=6, traininglength=120)
RMSFE = round(rbind(insample, out1ahead, out3head, out6ahead),3)
nnames=c(
  "VAR-BIC", "VAR-HQC", "AR-BIC", "AR-HQC",
  paste("VAR_",c(3,4,6),"(1)",sep=""),
  paste("AR_",c(3,4,6),"(1)",sep=""),
  "DNS-VAR", "DNS-AR"
)
colnames(RMSFE) = nnames
rownames(RMSFE) = c("in-sample 1-step", "1-ahead", "3-ahead", "6-ahead")
RMSFE
Sys.time()-start
write.table(RMSFE,file=paste("./table6-US.csv", sep=''),append=T,col.names=F,row.names=F)
