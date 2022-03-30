## ####################################################################
## ####################################################################
## Supplement for
## "Dynamic Factor Model for Functional Time Series: 
## Identification, Estimation, and Prediction"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Table 2.
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
library(parallel)
library(dffm)
## ##################################
## Cluster setup (ignore this when running on local machine)
## ##################################
if(is.na(strtoi(Sys.getenv(c("SLURM_NTASKS"))))){
  cl = makeCluster(detectCores()-1)
} else {
  ntasks <- strtoi(Sys.getenv(c("SLURM_NTASKS")))
  nslaves <- ntasks-1
  cl = makeCluster(nslaves, type="MPI")
}
input = as.numeric(commandArgs(trailingOnly = TRUE))
input
i = input[1]
type = input[2]
## ##################################
## Reproducible random number generator
## ##################################
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
snow::clusterSetupRNG(cl)
## ##################################
## Setup (for running on local machine)
## ##################################
if(is.na(i)) i=1 # i=1 for T=100, i=2 for T=200, i=3 for T=500
if(is.na(type)) type = 1 # type=1 for M1, type=2 for M2, etc.
MC=10 # number of Monte Carlo repetitions
## ##################################
T = c(100,200,500)[i]
trueK = c(3,2,2,1)
truep = c(1,2,4,4)
## ##################################
sim.criterion=function(T, type){
  fourierbasis = function(K, G = 21){
    # creates a list of the first K Fourier basis functions on [0,1] with gridsize G
    GRID = (0:(G-1))/(G-1)
    phi = list()
    phi[[1]] = rep(1, length(GRID))
    if(K > 1){
      for(i in 1:floor(K/2)){
        phi[[2*i]] = sqrt(2)*sin(2*i*pi*GRID)
        phi[[2*i+1]] = sqrt(2)*cos(2*i*pi*GRID)
      }
      phi[[K+1]] = NULL
    }
    return(phi)
  }
  gen.VARfactor = function(A, T = 100){
    # A is a quadratic coefficient matrix or a list of quadratic VAR coefficient matrices. The lag order is the length of the list
    # T is the sample size
    if( is.matrix(A) == TRUE ) A = list(A)
    if( is.list(A) != TRUE ) stop('input is not a list')
    if(length(unique(unlist(lapply(A, dim)))) != 1) stop('input matrices are either not quadratic or not compatible')
    p = length(A)
    K = dim(A[[1]])[1]
    f = matrix(0, nrow = K, ncol = T+p)
    eta = diag(1/(1:K))%*% matrix(rnorm(K*(T+p)), nrow = K, ncol = T+p)
    # eta = diag(seq(from=1, to=0, length.out = K)) %*% matrix(rnorm(K*(T+p)), nrow = K, ncol = T+p)
    for(t in (p+1):(p+T)){
      for(i in 1:p){
        f[,t] = f[,t] + A[[i]] %*% f[,t-i]
      }
      f[,t] = f[,t] + eta[,t]
    }
    factors = f[,-(1:p)]
    return(factors)
  }
  gen.FTS = function(A, T = 100, G = 21){
    # A is a quadratic coefficient matrix or a list of quadratic VAR coefficient matrices. The lag order is the length of the list
    # T is the sample size
    if( is.matrix(A) == TRUE ) A = list(A)
    if( is.list(A) != TRUE ) stop('input is not a list')
    if(length(unique(unlist(lapply(A, dim)))) != 1) stop('input matrices are either not quadratic or not compatible')
    factors = gen.VARfactor(A, T)
    K = dim(A[[1]])[1]
    basis = fourierbasis(K, G)
    Y = matrix(0, nrow = T, ncol = G)
    for(l in 1:K){
      Y = Y + t(factors[l,,drop=FALSE]) %*% t(basis[[l]])
    }
    colnames(Y) = round((0:(G-1))/(G-1),2)
    out = ts(Y, start = 2000, frequency = 12)
    return(out)
  }
  gen.data = function(type = 1, T = 100){
    if (type == 1){
      A1 = matrix(0, nrow = 10, ncol = 10)
      A1[1,1:3] = c(-0.05, -0.23, 0.76)
      A1[2,1:3] = c(0.80, -0.05, 0.04)
      A1[3,1:3] = c(0.04, 0.7, 0.23)
      data = gen.FTS(A1, T)
    } else if (type == 2) {
      A1 = matrix(0, nrow = 10, ncol = 10)
      A2 = A1
      A1[1,1:2] = c(0.8, -0.5)
      A1[2,1:2] = c(0.1, -0.5)
      A2[1,1:2] = c(-0.3, -0.3)
      A2[2,1:2] = c(-0.2, 0.3)
      data = gen.FTS(list(A1, A2), T)
    } else if (type == 3) {
      A1 = matrix(0, nrow = 10, ncol = 10)
      A2 = A1
      A3 = A1
      A4 = A1
      A1[1,1:2] = c(0.4, -0.2)
      A1[2,1:2] = c(0,0.3)
      A2[1,1:2] = c(-0.1, -0.1)
      A2[2,1:2] = c(0,-0.1)
      A3[1,1:2] = c(0.15,0.15)
      A3[2,1:2] = c(0,0.15)
      A4[1,1:2] = c(0.3, -0.4)
      A4[2,1:2] = c(0, 0.6)
      data = gen.FTS(list(A1,A2,A3,A4), T)
    } else if (type == 4) {
      A1 = matrix(0, nrow = 10, ncol = 10)
      A2 = A1
      A3 = A1
      A4 = A1
      A1[1,1] = 0.2
      A4[1,1] = 0.7
      data = gen.FTS(list(A1,A2,A3,A4), T)
    } else {
      stop('input is not valid')
    }
    return(data)
  }
  source("00_extrafunctions.R", print.eval = TRUE)
  data = gen.data(type = type, T)
  out=dffm.criterion(fpca.preprocess(data), K.max=8, p.max=8)
  return(out[[1]])
}
## ##################################
out1 = parSapply(cl, rep(T,MC), sim.criterion, type=type)
out.names = paste(rep(c("BIC", "HQC", "fFPE"),2), c(rep("-K",3), rep("-p",3)), sep="")
rownames(out1) = out.names
bias = rowMeans(out1) - c(rep(trueK[type],3), rep(truep[type],3))
rmse = sqrt(rowMeans(out1^2) - rowMeans(out1)^2 + bias^2)
bias
rmse
## ##################################
table.out = matrix(c(T,type,MC,bias,rmse), ncol=length(bias)+length(rmse)+3)
write.table(table.out,file=paste("./table2.csv", sep=''),append=T,col.names=F,row.names=F)
## ##################################
Sys.time()-start