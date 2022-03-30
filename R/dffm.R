#' dffm
#'
#' dffm is a function which performs a complete functional data analysis on a multivariate time series.
#'
#' @param fpcaobj A multivariate time series of class 'ts' or 'data.frame' or an object of the class 'FPCAobj'.
#' @param K A parameter which specifies the number of factors be used. If K is NULL, K will be 4. NULL is default.
#' @param p A parameter which specifies the number of lags will be used for the factor dynamics.
#' @param AR Will determine if factor dynamics will be based on vector autoregressive models or autoregressive models.
#' If AR is FALSE computations for factor dynamics will be based on vector autoregressive models. FALSE is default.
#'
#' @return
#' dffm will return an object of class 'dffm' with:
#' \item{K}{Used number of factors in the analysis.}
#' \item{p}{Used number of lags used for the factor dynamics.}
#' \item{fpcascores}{A matrix/ts-matrix of the computed factorscores.}
#' \item{factordynamics}{Used factor dynamics method.}
#' \item{VARcoefficients}{A matrix of the computed factor dynamics.}
#' \item{loadingfunctions.obsgrid}{A matrix containing the estimated loading functions, which are evaluated on the observationgrid.}
#' \item{meanfunction.obsgrid}{A vector containing the estimated mean function, which is evalauted on the observationgrid.}
#' \item{loadingfunctions.workgrid}{A matrix containing the estimated loading functions, which are evaluated on the observationgrid.}
#' \item{meanfunction.workgrid}{A vector containing the estimated mean function, which is evalauted on the workinggrid.}
#' \item{fittedcurve.obsgrid}{A matrix/ts-matrix of the fitted curves evaluated on the observationgrid.}
#' \item{fittedcurve.workgrid}{A matrix/ts-matrix of the fitted curves evaluated on the workinggrid.}
#' \item{observationgrid}{Used observationgrid in the analysis.}
#' \item{workinggrid}{Used workinggrid in the analysis.}
#' \item{FPCA}{Object of class 'FPCAobj' computed from the 'fpca.preprocess' function. For more details see 'fpca.preprocess' function.}
#'
#' @export
#' @references
#' * Kokoszka P. and Reimherr M. (2017), "Introduction to Functional Data Analysis". CRC Press.
#' * Otto S. and Salish N. (2022). "Dynamic Factor Model for Functional Time Series: Identification, Estimation, and Prediction". ...
#' * Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). Functional data analysis for sparse longitudinal
#' data. Journal of the American Statistical Association, 100:577-590.
#' @examples
#' # with object of class 'FPCAobj'
#' fed = load.fed()
#' fpca = fpca.preprocess(data = fed, method = "splines")
#' dffm(fpca)
#' # with raw data
#' fed = load.fed()
#' dffm(fed, K = 3, p = 2, AR = FALSE)
dffm = function(fpcaobj, K = NULL, p = NULL, AR = FALSE){
  # checking data
  if(class(fpcaobj)[1]!="FPCAobj"&class(fpcaobj)[1]=="mts"|class(fpcaobj)[1]!="FPCAobj"& class(fpcaobj)[1]=="data.frame"){
    fpcaobj = fpca.preprocess(fpcaobj)
  }
  if(class(fpcaobj)[1] != "FPCAobj" & class(fpcaobj)[1] != "mts" & class(fpcaobj)[1] != "data.frame"){
    stop("fpcaobj has to be an object of class FPCAobj or a multivariate time series of class 'ts' or 'data.frame'")
  }
  if(is.null(K)) K = 4
  if(is.null(p)) p = 1
  if(K > length(fpcaobj$eigenvalues)) stop(paste("K is not allowed to be greater then", length(fpcaobj$eigenvalues)))
  # VAR/AR matrix
  loadingfunctions.obsgrid = fpcaobj$eigenfunctions.obsgrid[,1:K]
  loadingfunctions.workgrid = fpcaobj$eigenfunctions.workgrid[,1:K]
  factors = fpcaobj$fpcascores[,1:K,drop=F]
  if(K == 1){
    factors = as.matrix(factors)
    if(is.ts(fpcaobj$raw.data)) factors = ts(factors, start = head(time(fpcaobj$raw.data),1), frequency = frequency(fpcaobj$raw.data))
  }
  if(AR==TRUE){
    VARcoefficients = matrix(0, nrow = K*p, ncol=K)
    for(i in 1:K){
      ARfit = lm(embed(factors[,i],p+1)[,1] ~ -1 + embed(factors[,i],p+1)[,-1])
      for(j in 1:p) VARcoefficients[i+K*(j-1),i] = ARfit$coefficients[j]
    }
  } else {
    VARcoefficients = lm(embed(factors, p+1)[,1:K] ~ -1 + embed(factors, p+1)[,-(1:K)])$coefficients
  }
  if(K > 1) colnames(VARcoefficients) = paste("F",1:K,sep="")
  if(K > 1) rownames(VARcoefficients) = paste("A",rep(1:p,each=K),"F",rep(1:K,p),sep="")
  # fitting curves
  temp.obsgrid = factors %*% t(loadingfunctions.obsgrid)
  fittedcurve.obsgrid = matrix(nrow = dim(temp.obsgrid)[1], ncol = dim(temp.obsgrid)[2])
  for (i in 1:dim(temp.obsgrid)[1]) {
    fittedcurve.obsgrid[i,] = temp.obsgrid[i,] + fpcaobj$meanfunction.obsgrid
  }
  temp.workgrid = factors %*% t(loadingfunctions.workgrid)
  fittedcurve.workgrid = matrix(nrow = dim(temp.workgrid)[1], ncol = dim(temp.workgrid)[2])
  for (i in 1:dim(temp.workgrid)[1]) {
    fittedcurve.workgrid[i,] = temp.workgrid[i,] + fpcaobj$meanfunction.workgrid
  }
  colnames(fittedcurve.obsgrid) = fpcaobj$observationgrid
  if(is.ts(fpcaobj$raw.data)) fittedcurve.obsgrid = ts(fittedcurve.obsgrid, start = head(time(fpcaobj$raw.data),1), frequency = frequency(fpcaobj$raw.data))
  colnames(fittedcurve.workgrid) = fpcaobj$workinggrid
  if(is.ts(fpcaobj$raw.data)) fittedcurve.workgrid = ts(fittedcurve.workgrid, start = head(time(fpcaobj$raw.data),1), frequency = frequency(fpcaobj$raw.data))
  out=list(
    "K" = K,
    "p" = p,
    "factors" = factors,
    "factordynamics" = c("VAR","AR")[AR+1],
    "VARcoefficients" = VARcoefficients,
    "loadingfunctions.obsgrid" = loadingfunctions.obsgrid,
    "meanfunction.obsgrid" = fpcaobj$meanfunction.obsgrid,
    "loadingfunctions.workgrid" = loadingfunctions.workgrid,
    "meanfunction.workgrid" = fpcaobj$meanfunction.workgrid,
    "fittedcurve.obsgrid" = fittedcurve.obsgrid,
    "fittedcurve.workgrid" = fittedcurve.workgrid,
    "observationgrid" = fpcaobj$observationgrid,
    "workinggrid" = fpcaobj$workinggrid,
    "FPCA" = fpcaobj
  )
  class(out) = "dffm"
  return(out)
}

#' print.dffm
#'
#' S3 print function for dffm objects.
#'
#' @param x An object of class "dffm".
#'
#' @return
#' @export
#'
#' @examples
print.dffm = function(x){
  print(list(
    "K" = x$K,
    "p" = x$p,
    "factors" = paste(paste("class:", class(x$factors)[1],";",dim(x$factors)[1],
                            "rows,", dim(x$factors)[2], "cols")),
    "factordynamics" = x$factordynamics,
    "VARcoefficients" = paste(paste("class:", class(x$VARcoefficients)[1],";",dim(x$VARcoefficients)[1],
                                    "rows,", dim(x$VARcoefficients)[2], "cols")),
    "loadingfunctions.obsgrid" = paste(paste("class:", class(x$loadingfunctions.obsgrid)[1],";",dim(x$loadingfunctions.obsgrid)[1],
                                             "rows,", dim(x$loadingfunctions.obsgrid)[2], "cols")),
    "loadingfunctions.obsgrid" = paste(paste("class:", class(x$loadingfunctions.obsgrid)[1],";",dim(x$loadingfunctions.obsgrid)[1],
                                             "rows,", dim(x$loadingfunctions.obsgrid)[2], "cols")),
    "loadingfunctions.workgrid" = paste(paste("class:", class(x$loadingfunctions.workgrid)[1],";",dim(x$loadingfunctions.workgrid)[1],
                                             "rows,", dim(x$loadingfunctions.workgrid)[2], "cols")),
    "meanfunction.workgrid" = paste(paste("class:", class(x$meanfunction.workgrid)[1],";",dim(x$meanfunction.workgrid)[1],
                                             "rows,", dim(x$meanfunction.workgrid)[2], "cols")),
    "fittedcurve.obsgrid" = paste(paste("class:", class(x$fittedcurve.obsgrid)[1],";",dim(x$fittedcurve.obsgrid)[1],
                                          "rows,", dim(x$fittedcurve.obsgrid)[2], "cols")),
    "fittedcurve.workgrid" = paste(paste("class:", class(x$fittedcurve.workgrid)[1],";",dim(x$fittedcurve.workgrid)[1],
                                          "rows,", dim(x$fittedcurve.workgrid)[2], "cols")),
    "observationgrid" = x$observationgrid,
    "workinggrid" = x$workinggrid,
    "FPCA" = paste("class:", class(x$FPCA))
  ), quote = FALSE)
}

#' predict.dffm
#'
#' S3 predict function for dffm objects.
#'
#' @param x An object of class "dffm".
#' @param h Time ahead forecasts.
#'
#' @return
#' @export
#'
#' @examples
predict.dffm = function(x, h = NULL){
  if(is.null(h)) h = 1
  output = dffm.forecast(dffmobj = x, h = h)
  return(output)
}

#' plot.dffm
#'
#' S3 predict function for dffm objects.
#'
#' @param x An object of class "dffm".
#' @param time Observation to be plotted.
#' @param date.on Logical; if TRUE, date will be shown underneath plot. If FALSE, no date will be shown. FALSE is default.
#' @param workgrid Logical; If TRUE, workgrid will be used. If FALSE observationgrid will be used. TRUE is default.
#' @param ... Plot adjustment parameters left open.
#'
#' @return
#' @export
#'
#' @examples
plot.dffm = function(x, time = NULL, date.on = FALSE, workgrid = TRUE, ...){
  if(is.null(time)) time = 1
  if(date.on == TRUE & class(x$fittedcurve.obsgrid)[1] != "mts") stop("'dffm' object has no time series 'fitted.data'")
  if(workgrid == TRUE){
    grid = as.numeric(x$workinggrid)
    y = x$fittedcurve.workgrid[time,]
  }else{
    grid = as.numeric(x$observationgrid)
    y = x$fittedcurve.obsgrid[time,]
  }
  if(class(x$fittedcurve.obsgrid)[1] == "mts" & date.on == TRUE){
    plot(x = grid, y = y, type = "l",
         sub = list(paste("from", zoo:: as.yearmon(time(x$fittedcurve.obsgrid))[time]), cex =.9), ...)
  }
  if(class(x$fittedcurve.obsgrid)[1] == "mts" & date.on == FALSE){
    plot(x = grid, y = y, type = "l", ...)
  }
  if(class(x$fittedcurve.obsgrid)[1] != "mts"){
    plot(x = grid, y = y, type = "l", ...)
  }
}









