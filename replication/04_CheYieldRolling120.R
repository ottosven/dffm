rm(list=ls())

## run time approximately 4 hours

library(dffm)
library(glmnet)

data = readRDS("./data/LW360.rds")

get.results = function(index, data, rolling = TRUE, h=1, wsize = 120){
  start = 1
  if(rolling) start = index
  end = index + wsize - 1
  training = data[start:end,]
  testing = data[end+h,]
  
  fdaobj = fda.preprocess(training)
  cumacobj = fts.cumAC(fdaobj)
  IC = fts.criterion(fdaobj, K.max = 10, p.max = 10)$IC.min
  
  #### ##########################################################

  shrinkageforecast = function(obj, K, p, h=1){
    factors = obj$scores.centered[,1:K, drop=F]
    ts_data = factors
    prepare_lagged_data <- function(ts_data, p) {
      embedded_data <- embed(ts_data, p + 1)  # Creates lagged matrix
      X <- embedded_data[, -c(1:ncol(ts_data)), drop = FALSE]  # Past p lags as features
      Y <- embedded_data[, 1:ncol(ts_data), drop = FALSE]  # Current values as targets
      list(X = as.matrix(X), Y = as.matrix(Y))
    }
    data_lagged <- prepare_lagged_data(ts_data, p = p)
    X_train <- data_lagged$X
    Y_train <- data_lagged$Y
    X_test = matrix(tail(embed(factors, p), 1), nrow=1)
    
    lassopred = function(...){
      cv_glmnet_models <- lapply(1:K, function(i) {
        cv.glmnet(X_train, Y_train[, i], alpha = 1, nfolds = 10)  # 10-fold CV
      })
      # Extract best lambda for each model
      best_lambdas <- sapply(cv_glmnet_models, function(model) model$lambda.min)
      glmnet_models <- lapply(1:K, function(i) {
        glmnet(X_train, Y_train[, i], alpha = 1, lambda = best_lambdas[i])
      })
      predictors = X_test
      factors.predict <- sapply(glmnet_models, function(model) {
        predict(model, newx = predictors)
      })
      if(h > 1){
        for(i in 2:h){
          if(p > 0) {
            predictors = matrix(c(factors.predict, predictors[1:(K*(max(p-1, 0)))]), nrow=1)
          } else {
            predictors = matrix(factors.predict, nrow=1)
          }
          factors.predict <- sapply(glmnet_models, function(model) {
            predict(model, newx = predictors)
          })
        }
      }
      return(factors.predict)
    }
    
    ridgepred = function(...){
      cv_glmnet_models <- lapply(1:K, function(i) {
        cv.glmnet(X_train, Y_train[, i], alpha = 0, nfolds = 10)  # 10-fold CV
      })
      # Extract best lambda for each model
      best_lambdas <- sapply(cv_glmnet_models, function(model) model$lambda.min)
      glmnet_models <- lapply(1:K, function(i) {
        glmnet(X_train, Y_train[, i], alpha = 0, lambda = best_lambdas[i])
      })
      predictors = X_test
      factors.predict <- sapply(glmnet_models, function(model) {
        predict(model, newx = predictors)
      })
      if(h > 1){
        for(i in 2:h){
          if(p > 0) {
            # Fix: Ensure proper matrix structure is maintained
            predictors = matrix(c(factors.predict, predictors[1:(K*(max(p-1, 0)))]), nrow=1)
          } else {
            predictors = matrix(factors.predict, nrow=1)
          }
          factors.predict <- sapply(glmnet_models, function(model) {
            predict(model, newx = predictors)
          })
        }
      }
      return(factors.predict)
    }
    loadings = obj$eigenfunctions[,1:K]
    factors.lasso = lassopred()
    lasso = factors.lasso %*% t(loadings) + obj$meanfunction
    factors.ridge = ridgepred()
    ridge = factors.ridge %*% t(loadings) + obj$meanfunction
    return(t(rbind(lasso, ridge)))
  }
  
  RW = data[end,,drop=F]
  DNSRW = t(rbind(
    RW,
    DNS.forecast(fdaobj, p=1, AR = FALSE, h=h),
    DNS.forecast(fdaobj, p=2, AR = FALSE, h=h)
  ))
  colnames(DNSRW) = c("RW", "DNS.VAR1", "DNS.VAR2")

  BIC = fts.VARforecast(cumacobj, K = IC[1,1], p = IC[1,2], AR = FALSE, h=h)
  HQC = fts.VARforecast(cumacobj, K = IC[2,1], p = IC[2,2], AR = FALSE, h=h)
  FPC.fFPE = fts.VARforecast(fdaobj, K = IC[3,1], p = IC[3,2], AR = FALSE, h=h)
  
  ICbased = cbind(t(rbind(BIC, HQC, FPC.fFPE)))
  colnames(ICbased) = c("BIC.ols", "HQC.ols", "FPC.ols")
  
  BIC.shrink = shrinkageforecast(cumacobj, K=IC[1,1], p=IC[1,2])
  HQC.shrink = shrinkageforecast(cumacobj, K=IC[2,1], p=IC[2,2])
  FPC.fFPE.shrink = shrinkageforecast(fdaobj, K=IC[3,1], p=IC[3,2])
  
  ICbased.shrink = cbind(BIC.shrink, HQC.shrink, FPC.fFPE.shrink)
  
  colnames(BIC.shrink) = c("BIC.lasso", "BIC.ridge")
  colnames(HQC.shrink) = c("HQC.lasso", "HQC.ridge")
  colnames(FPC.fFPE.shrink) = c("FPC.lasso", "FPC.ridge")
  
  VARforecasts = function(K, obj, AR=FALSE, p=1){
    fts.VARforecast(obj, K = K, p = p, AR = AR, h=h)
  }
  VAR1.allK = cbind(cumacobj$meanfunction, 
                  sapply(1:15, VARforecasts, obj = cumacobj, AR=FALSE, p=1))
  colnames(VAR1.allK) = paste0("VAR1.K",0:15)
  FPC.allK = cbind(cumacobj$meanfunction, 
                    sapply(1:15, VARforecasts, obj = fdaobj, AR=FALSE, p=1))
  colnames(FPC.allK) = paste0("FPC.K",0:15)
  
  #### ##########################################################
  
  allforecasts = cbind(DNSRW, ICbased, ICbased.shrink, VAR1.allK, FPC.allK)
  
  forecasterrors = (allforecasts - testing)^2
  return(forecasterrors)  
}

## ###############################################################################
## ###############################################################################
start = Sys.time()
wsize = 120
## ###############################################################################

mypath = "./CheResults/YieldRolling120-"


h=1
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "1step.rds"))
Sys.time() - start

h=3
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "3step.rds"))

h=6
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "6step.rds"))

h=12
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "12step.rds"))


Sys.time() - start
