## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Figure 1.
## ####################################################################
## ####################################################################
library(dffm)

MSEmatrix = as.matrix(read.csv("./replication/MSEmatrix-example.csv"))
IC.example = list()
IC.example[["MSE.matrix"]] = MSEmatrix
## ##################################
## Plot Figure 1
## ##################################
pdf("figure1.pdf", width=18, height=16, pointsize = 30)
dffm.MSEplot(IC.example, kmax=8, pmax=6)
dev.off()
