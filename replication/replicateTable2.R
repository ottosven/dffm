## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Table 2.
## ####################################################################
## ####################################################################
rm(list=ls())
library(dffm)
library(tidyverse)

get.RMSEdata = function(type, n){
  file_list = list.files(
    path = "./replication/simulationresults",
    pattern = paste("^ICsim-",type,"-",n,sep=""))
  data_frames <- list()
  # Loop over the file list, read each file and store it in the list
  for (file_name in file_list) {
    data_frame = read.table(paste("./replication/simulationresults/",file_name,sep=""))
    # data_frame <- read.csv(file_name, stringsAsFactors = FALSE) # adjust read.csv parameters as needed
    data_frames[[length(data_frames) + 1]] <- data_frame
  }
  # Combine all data frames into one large data frame
  combined_data_frame <- do.call(rbind, data_frames)
  # Data on information criteria
  ICdata = combined_data_frame |> select(V4,V5,V6,V7)
  if(type %in% c(1,2,3)){
    trueparam = c(1,1,2,2)
  } else if(type == 4){
    trueparam = c(2,2,3,3)
  } else {
    trueparam = c(3,3,1,1)
  }

  est.error = sweep(ICdata, MARGIN = 2, STATS = trueparam, FUN = function(x, y) y - x)
  falsesel = colMeans(ifelse(est.error == 0, FALSE, TRUE))
  bias = colMeans(est.error)
  RMSE = sqrt(colMeans(ICdata^2) - colMeans(ICdata)^2 + bias^2)
  # Data on chi-hat
  chihat.errornorm = combined_data_frame |> select(V8) |> colMeans()
  return(c(chihat.errornorm, RMSE, bias, falsesel))
}

## ####################################################
## Evaluate
## ####################################################

results.table = rbind(
  get.RMSEdata(1,100),
  get.RMSEdata(1,200),
  get.RMSEdata(1,500),
  get.RMSEdata(2,100),
  get.RMSEdata(2,200),
  get.RMSEdata(2,500),
  get.RMSEdata(3,100),
  get.RMSEdata(3,200),
  get.RMSEdata(3,500),
  get.RMSEdata(4,100),
  get.RMSEdata(4,200),
  get.RMSEdata(4,500),
  get.RMSEdata(5,100),
  get.RMSEdata(5,200),
  get.RMSEdata(5,500)
)

colnames(results.table) = c("chihat.meanerr",
                            "RMSE.Kbic", "RMSE.Khqc", "RMSE.pbic", "RMSE.phqc",
                            "bias.Kbic", "bias.Khqc", "bias.pbic", "bias.phqc",
                            "false.Kbic", "false.Khqc", "false.pbic", "false.phqc"
)
rownames(results.table) = c(
  paste0("M1T",c(100,200,500)), paste0("M2T",c(100,200,500)),
  paste0("M3T",c(100,200,500)),
  paste0("M4T",c(100,200,500)), paste0("M5T",c(100,200,500))
)

results = round(results.table,2)
results
write.table(results,file=paste("./table2.csv", sep=''),append=T,col.names=F,row.names=F)
