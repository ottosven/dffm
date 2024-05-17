## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Figure 1.
## ####################################################################
## ####################################################################
rm(list=ls())
library(dffm)
library(tidyverse)

get.psihatplot = function(type, n){
  file_list = list.files(
    path = "./replication/simulationresults",
    pattern = paste("^PSIsim-",type,"-",n,sep=""))
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
  PSIdata = combined_data_frame |> select(-c(1:3))

  G = dim(PSIdata)[2]
  grid = seq(0,1,length.out = G)
  quantiles_5_95 = apply(PSIdata, 2, quantile, probs = c(0.05, 0.95))
  if(type == 1){
    plot(grid, quantiles_5_95[1,], lty = 2, type="l", ylim = c(0.1,1.9), xlab="", ylab="")
  } else {
    plot(grid, quantiles_5_95[1,], lty = 2, type="l", ylim = c(-2.4,2.4), xlab="", ylab="")
  }
  lines(grid, quantiles_5_95[2,], lty = 2)
  polygon(c(grid, rev(grid)), c(quantiles_5_95[1,], rev(quantiles_5_95[2,])), col = "lightgrey", border = NA, density = 10)
  lines(grid, colMeans(PSIdata))
  if(type == 1){
    lines(grid, rep(1,G), lty = 3, lwd = 6)
  } else if(type == 2){
    lines(grid, sqrt(2)*sin(2*pi*grid), lty = 3, lwd = 6)
  } else {
    lines(grid, sqrt(2)*sin(4*pi*grid), lty = 3, lwd = 6)
  }
}


pdf("./figure2.pdf", width=30, height=15, pointsize = 30)
par(mfrow = c(3,3), cex=1, mai = c(1, 1, 1, 1))
get.psihatplot(1,100)
get.psihatplot(2,100)
get.psihatplot(3,100)
get.psihatplot(1,200)
get.psihatplot(2,200)
get.psihatplot(3,200)
get.psihatplot(1,500)
get.psihatplot(2,500)
get.psihatplot(3,500)
dev.off()
