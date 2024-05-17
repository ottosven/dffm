## ####################################################################
## ####################################################################
## Supplement for
## "Approximate Factor Models for Functional Time Series"
## by Sven Otto and Nazarii Salish.
## This R-script allows to reproduce Figure 5.
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
## For each time series use the following parameters:
## Start:01.01.1995, End:30.06.2022, Interval:1MO
## 4) Save the data for each country into a mts-object with monthly frequency. The column names are the times to maturity in months.
## Note that some values are missing in the original dataset. Missing vallues in the mts-object are set to "NA".
## 5) Save the list of the six mts-objects as "G7data.Rdata". The names of the list elements are c("CA", "FR", "DE", "IT", "JP", "GB")
## ####################################################################
## ####################################################################
library(dffm)
## ##################################
## Data
## ##################################
load("./not-for-public-access/eikonyields/G7data.RData")
## ##################################
## Information Criteria and loading functions
## ##################################
get.workinggrid = function(data, gridsize){
  obsgrid = as.numeric(colnames(data))
  seq(obsgrid[1], obsgrid[length(obsgrid)], gridsize)
}
##
figure.loadingfunctions = function(data, K = 4, signpattern = NULL, main = NULL, ylim = NULL){
  if(is.null(signpattern)) signpattern=rep(1,K)
  if(is.null(main)) main=paste("Loadings", deparse(substitute(data)), "1995-2022")
  FDAdata = fda.preprocess(data, workinggrid = get.workinggrid(data, 0.5))
  loadings = FDAdata$eigenfunctions[,1:K]
  loadings = loadings%*%diag(signpattern)
  if(is.null(ylim)) ylim = range(FPCloadings)
  plot(FDAdata$workinggrid, loadings[,1], ylim = ylim, type='l', lty=1, col=1, lwd=3,
       xlab = 'Maturity (months)', ylab = '', main=main,  cex.lab=1.3, cex.axis=1.3, cex.main=1.9, yaxt="n")
  axis(2, at = seq(-0.1, 0.2, by = 0.1), cex.axis=1.3)
  for(i in 2:K){
    lines(FDAdata$workinggrid, loadings[,i], lty=i, col=1, lwd = 4-i+3)
  }
  legendnames = c(expression(hat(psi)[1]),expression(hat(psi)[2]), expression(hat(psi)[3]), expression(hat(psi)[4]), expression(hat(psi)[5]), expression(hat(psi)[6]))
  legend(x="top",
         legend = legendnames[1:K],
         col = 1,
         lty = 1:K,
         lwd = c(3,5,4,3,2,1)[1:K],
         pt.cex = 2,
         cex = 1,
         text.col = "black",
         seg.len = 2.5,
         horiz = TRUE,
         bty = "n"
  )
}
## ##################################
## Plot loading functions
## ##################################
pdf("figure5.pdf", width=30, height=18, pointsize = 30)
par(mfrow=c(2,3))
figure.loadingfunctions(G7data$CA, 4, signpattern = c(1,1,-1,-1), ylim = c(-0.12, 0.188), main = "Canada")
figure.loadingfunctions(G7data$FR, 4, signpattern = c(1,1,-1,-1), ylim = c(-0.12, 0.188), main = "France")
figure.loadingfunctions(G7data$DE, 4, signpattern = c(1,1,-1,-1), ylim = c(-0.12, 0.188), main = "Germany")
figure.loadingfunctions(G7data$IT, 4, signpattern = c(1,1,-1,-1), ylim = c(-0.12, 0.188), main = "Italy")
figure.loadingfunctions(G7data$JP, 4, signpattern = c(1,1,-1,-1), ylim = c(-0.12, 0.188), main = "Japan")
figure.loadingfunctions(G7data$GB, 4, signpattern = c(1,1,-1,-1), ylim = c(-0.12, 0.188), main = "United Kingdom")
dev.off()
