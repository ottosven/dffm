## ####################################################################
## ####################################################################
## Supplement for
## "Dynamic Factor Model for Functional Time Series:
## Identification, Estimation, and Prediction"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Table 4.
## ####################################################################
## ####################################################################
#####################
## INFO: Code runs about 20-30 minutes
#####################
library(dffm)
## ##################################
## Data
## ##################################
fed = load.fed()
JKV = load.JKV()
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
insample.prederror = function(data){
  prederrors = list()
  for(j in 1:12) prederrors[[j]] = matrix(nrow = dim(data)[1]-1, ncol = dim(data)[2])
  data.fpca = fpca.preprocess(data)
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
  for(j in 1:12) rmse[j] = sqrt(mean(prederrors[[j]]^2))
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
    data.fpca = fpca.preprocess(traindata)
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
  for(j in 1:12) rmse[j] = sqrt(mean(prederrors[[j]]^2))
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
## Results for JKV
## ##################################
insample = insample.prederror(JKV)
out1ahead = outofsample.prederror(JKV, h=1, traininglength=120)
out3head = outofsample.prederror(JKV, h=3, traininglength=120)
out6ahead = outofsample.prederror(JKV, h=6, traininglength=120)
out12ahead = outofsample.prederror(JKV, h=12, traininglength=120)
RMSFE.JKV = round(rbind(insample, out1ahead, out3head, out6ahead, out12ahead),3)
nnames=c(
  "VAR-BIC", "VAR-HQC", "AR-BIC", "AR-HQC",
  paste("VAR_",c(3,4,6),"(1)",sep=""),
  paste("AR_",c(3,4,6),"(1)",sep=""),
  "DNS-VAR", "DNS-AR"
)
colnames(RMSFE.JKV) = nnames
rownames(RMSFE.JKV) = c("in-sample 1-step", "1-ahead", "3-ahead", "6-ahead", "12-ahead")
RMSFE.JKV
## ##################################
## Results for fed
## ##################################
insample = insample.prederror(fed)
out1ahead = outofsample.prederror(fed, h=1, traininglength=120)
out3head = outofsample.prederror(fed, h=3, traininglength=120)
out6ahead = outofsample.prederror(fed, h=6, traininglength=120)
out12ahead = outofsample.prederror(fed, h=12, traininglength=120)
RMSFE.fed = round(rbind(insample, out1ahead, out3head, out6ahead, out12ahead),3)
nnames=c(
  "VAR-BIC", "VAR-HQC", "AR-BIC", "AR-HQC",
  paste("VAR_",c(3,4,6),"(1)",sep=""),
  paste("AR_",c(3,4,6),"(1)",sep=""),
  "DNS-VAR", "DNS-AR"
)
colnames(RMSFE.fed) = nnames
rownames(RMSFE.fed) = c("in-sample 1-step", "1-ahead", "3-ahead", "6-ahead", "12-ahead")
RMSFE.fed
