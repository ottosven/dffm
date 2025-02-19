library(dffm)

# approximately 30 minutes runtime

fdaobj = readRDS("./data/LW360-fda.rds")
fts.criterion(fdaobj, K.max = 15, p.max = 4)$IC

start = Sys.time()  
rollingK = function(data, rollingindices, width = 120){
  getIC = function(index, data){
    start = index
    end = index + width
    training = data[start:end,]
    fdaobj = fda.preprocess(training, workinggrid = seq(1, 360, by=1))
    crit = fts.criterion(fdaobj, K.max = 15, p.max = 4)
    c(crit$IC.min[1:2,1], crit$IC.min[1:2,2])
  }
  result = sapply(rollingindices, getIC, data = data)
  rownames(result) = c("BIC-K", "HQC-K", "BIC-p", "HQC-p")
  result
}
LW360 = readRDS("./data/LW360.rds")
data = LW360
rollingindices = 1:(dim(data)[1]-120)
out = rollingK(data = data, rollingindices = rollingindices, width = 120)
allK = ts(t(out), start = time(data)[1], frequency = 12)
allK
saveRDS(allK, "./yield-results/rollingK-w120.rds")
Sys.time() - start  