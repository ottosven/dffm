library(dffm)

start = Sys.time()

## run time approximately 4 minutes

data = readRDS("./data/mort.rds")

get.results = function(index, data, rolling = TRUE, h=1, wsize = 120){
  start = 1
  if(rolling) start = index
  end = index + wsize - 1
  training = data[start:end,]
  testing = data[end+h,]
  
  fdaobj = fda.preprocess(training)
  cumacobj = fts.cumAC(fdaobj)
  
  IC = fts.criterion(fdaobj, K.max = 10, p.max = 4)$IC.min
  BIC.VAR.AC = fts.VARforecast(cumacobj, K = IC[1,1], p = IC[1,2], AR = FALSE, h=h)
  HQC.VAR.AC = fts.VARforecast(cumacobj, K = IC[2,1], p = IC[2,2], AR = FALSE, h=h)
  fFPE.VAR.AC = fts.VARforecast(cumacobj, K = IC[3,1], p = IC[3,2], AR = FALSE, h=h)
  ICbased.AC = cbind(t(rbind(BIC.VAR.AC, HQC.VAR.AC, fFPE.VAR.AC)))
  colnames(ICbased.AC) = c("BIC.VAR.AC", "HQC.VAR.AC", "fFPE.VAR.AC")
  
  BIC.VAR.PC = fts.VARforecast(fdaobj, K = IC[1,1], p = IC[1,2], AR = FALSE, h=h)
  HQC.VAR.PC = fts.VARforecast(fdaobj, K = IC[2,1], p = IC[2,2], AR = FALSE, h=h)
  fFPE.VAR.PC = fts.VARforecast(fdaobj, K = IC[3,1], p = IC[3,2], AR = FALSE, h=h)
  ICbased.PC = cbind(t(rbind(BIC.VAR.PC, HQC.VAR.PC, fFPE.VAR.PC)))
  colnames(ICbased.PC) = c("BIC.VAR.PC", "HQC.VAR.PC", "fFPE.VAR.PC")
  
  RW = matrix(data[end,], ncol = 1)
  colnames(RW) = "RW"
  MF = matrix(colMeans(data), ncol = 1)
  colnames(MF) = "MF"

  ARforecasts = function(K, obj, AR=TRUE, p=1){
    fts.VARforecast(obj, K = K, p = p, AR = AR, h=h)
  }
  VAR1.AC = sapply(1:19, ARforecasts, obj = cumacobj, AR=FALSE, p=1)
  colnames(VAR1.AC) = paste0("VAR1-K",1:19,".AC")
  VAR2.AC = sapply(1:19, ARforecasts, obj = cumacobj, AR=FALSE, p=2)
  colnames(VAR2.AC) = paste0("VAR2-K",1:19,".AC")
  VAR1.PC = sapply(1:19, ARforecasts, obj = fdaobj, AR=FALSE, p=1)
  colnames(VAR1.PC) = paste0("VAR1-K",1:19,".PC")
  VAR2.PC = sapply(1:19, ARforecasts, obj = fdaobj, AR=FALSE, p=2)
  colnames(VAR2.PC) = paste0("VAR2-K",1:19,".PC")
  
  allpredictions = cbind(RW, MF, VAR1.AC, VAR2.AC, VAR1.PC, VAR2.PC, ICbased.AC, ICbased.PC)
  
  forecasterrors = (allpredictions - testing)^2
  return(forecasterrors)  
}

## ###############################################################################
## ###############################################################################
wsize = 50
mypath = "./mortality-results/"
## ###############################################################################

h=1
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "1step-rolling.rds"))

h=2
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "2step-rolling.rds"))

h=3
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "3step-rolling.rds"))

h=4
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "4step-rolling.rds"))

h=5
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "5step-rolling.rds"))

h=6
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "6step-rolling.rds"))

h=7
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "7step-rolling.rds"))

h=8
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "8step-rolling.rds"))

h=9
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "9step-rolling.rds"))

h=10
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "10step-rolling.rds"))

h=11
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "11step-rolling.rds"))

h=12
rollingindices = 1:(dim(data)[1]-wsize-h)
out = lapply(rollingindices, get.results, data = data, h=h, wsize=wsize, rolling=TRUE)
MSFEs = Reduce("+", out)/length(rollingindices)
saveRDS(MSFEs, paste0(mypath, "12step-rolling.rds"))

Sys.time() - start

