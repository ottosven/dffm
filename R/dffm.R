#' @import stats
#' @import fda
#' @import xts
#' @import tsbox
#' @import fdapace
#' @import vars
#' @import zoo
NULL

#' dffm
#'
#' dffm is a function which accompanies the fully functional factor model method from Otto and Salish (2021).
#' Furthermore does dffm perform the the pace method for sparse data from Yao et al. (2005) if needed.
#'
#' @param X A multivariate time series of class 'ts' or 'data.frame', an object of the class 'fd' from the package 'fda',
#' or an object of the class 'FPCA' from the package 'fdapace'.
#' @param pace Preprocessing method to be used. If TRUE the PACE method from the package 'fdapace' is applied (see Yao et al., 2005).
#' In case of irregularly observed functional data please set to TRUE.
#' If FALSE smoothing splines from the package 'fda' are used (see Kokoszka and Reimherr, 2017). FALSE is default.
#' @param criterion Automatic selection of the number of factors K. If TRUE the consistent information criterion of Otto and Salish (2021)
#' is used. If FALSE the value specified in 'K' is used. FALSE is default. For more details see 'dffm.criterion'.
#' @param observationgrid An optional list, sequence or numeric vector which contains a predetermined observationgrid for the
#' observations. If NULL the column names or input observationgrid from the 'fd' or 'FCPA' object will be used.
#' NULL is default.
#' @param user.gridsize An optional parameter specifying the size of the workinggrid. If user.gridsize is NULL, user.gridsize will be the highest
#' number from the column names or the specified observationgrid. NULL is default.
#' @param K A parameter which specifies the number of factors be used. If K is NULL, K will be 3. NULL is default.
#'
#' @return
#' dffm will return an object of class 'dffm' with:
#' \item{basis}{A list containing the call of the function; input observationgrid; used workinggrid; if cirterion is TRUE, the
#' optimal number of lags for prediction p; the number of used components K and input data if it was a data frame or ts object.}
#' \item{eigenvalues}{A vector containing the first K eigenvalues of the sample covariance operator.}
#' \item{loadingfunctions}{A matrix containing the estimated loading functions, which are evaluated on the workinggrid.}
#' \item{meanfunction}{A vector containing the estimated mean function, which is evalauted on the workinggrid.}
#' \item{factorscores}{A matrix/ts-matrix of the computed factorscores.}
#' \item{fitted.data}{A matrix/ts-matrix of the fitted curves evaluated on the workinggrid.}
#'
#' @export
#' @references
#' * Kokoszka P. and Reimherr M. (2017), "Introduction to Functional Data Analysis". CRC Press.
#' * Otto S. and Salish N. (2021). "Dynamic Factor Model for Functional Time Series: Identification, Estimation, and Prediction". ...
#' * Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). Functional data analysis for sparse longitudinal
#' data. Journal of the American Statistical Association, 100:577-590.
#' @seealso
#' * summary: summary of the dffm object.
#' * print: print of the dffm object.
#' * predict: predict of the dffm object. For more details, see dffm.predict.
#' * plot: plot of the dffm object at specified time with or without date.
#' * fitted: output of the fitted values of the dffm object.
#' @examples
#' # standard case
#' dl = load.dieboldli()
#' x = dffm(X = dl, pace = FALSE, criterion = FALSE, observationgrid = NULL, user.gridsize = NULL, K = 3)
#' x
#' summary(x)
#' plot(x, xlab = "maturities", ylab = "yields", date.on = TRUE, time = 20)
#' dffm.criterion(x = x) # optimal K = 3 and p = 1
#' p = predict(x = x, AR = FALSE, criterion = FALSE, h = 10, p = 1, K = 3)
#' plot(y = p$predicted.values[3,], x = 1:120, type = "l", xlab = "maturities", ylab = "yields", main = "predicted yields mar 2001")
#' # using dffm.preprocessing
#' pp = dffm.preprocessing(data = dl, pace = FALSE, lambda = 10^{-6})
#' x = dffm(X = pp, criterion = TRUE)
#' x$basis$p.opt  # optimal p = 1
#' x$basis$K  # optimal K = 3
#' predict(x, p = 1, K = 3, h = 3)
dffm = function(X, pace = FALSE, criterion = FALSE, observationgrid = NULL, user.gridsize = NULL, K = NULL){
  ## checking input information
  if(missing(X)) stop("please provide valid data or object. For more details see 'X' in the help file")
  if(class(X)[1] != "xts" & class(X)[1] != "mts" & class(X)[1] != "fd" & class(X)[1] != "FPCA" &
     class(X)[1] != "data.frame") stop("please provide valid data or object. For more details see 'X' in the help file")
  if(criterion == TRUE & !is.null(K))stop("please specify either K or criterion")
  if(class(X)[1] != "fd" & class(X)[1] != "FPCA"){
    if(sum(is.na(X[,1])) > 0 | sum(is.na(X[,dim(X)[2]])) > 0)stop("first and last column of X has to be without NA")}
  if(!is.null(K)){
    if(class(X)[1] == "mts" | class(X)[1] == "xts"){
      if(K > dim(X)[2]) stop("K is larger than the available number of components")
    }
    if(class(X)[1] == "fd"){
      if(K > X$basis$nbasis) stop("K is larger than the available number of components")
    }
    if(class(X)[1] == "FPCA"){
      if(K > length(X$lambda)) stop("K is larger than the available number of components")
    }
  }
  p = NULL
  X.name = substitute(X)
  obse.name = observationgrid
  grid.name = user.gridsize
  K.name = K
  criterion.name = criterion
  ## criterion
  if(criterion == TRUE){
    ## fd object
    if(class(X)[1] == "fd"){
      observationgrid = as.numeric(X$fdnames$time)
      if(is.null(user.gridsize)) user.gridsize = tail(observationgrid,1)
      if(sum(is.na(X$coefs))>0){
        stop("can not compute dffm because of NA in the 'fd' object. The pace method will not omit NA in the FPCA object")}
      g = seq(from = head(observationgrid,1), to = tail(observationgrid,1), by = tail(observationgrid,1)/user.gridsize)
      grid = g / tail(g,1)
      YieldPCA = pca.fd(X, nharm = X$basis$nbasis, harmfdPar = fdPar(X), centerfns = TRUE)
      Fhat <- YieldPCA$scores
      eigenvalues = YieldPCA$values[1:length(observationgrid)]
      observationgrid = obse.name
      user.gridsize = grid.name
    }
    ## FPCA object
    if(class(X)[1] == "FPCA"){
      observationgrid = as.numeric(X$fdnames$time)
      user.gridsize = tail(observationgrid,1)-head(observationgrid,1)+1
      Fhat = X$xiEst
      eigenvalues = X$lambda
      observationgrid = obse.name
      user.gridsize = grid.name
    }
    ## ts
    #fd
    if(class(X)[1] == "xts" & pace == FALSE | class(X)[1] == "mts" & pace == FALSE){
      if(is.null(observationgrid))  observationgrid = as.numeric(colnames(X))
      if(length(observationgrid) != ncol(X)) stop('please specify a valid observationgrid')
      if(is.na(as.numeric(observationgrid))[1])stop("column names have to be numbers or you have to specify a observationgrid")
      observationgrid = as.numeric(observationgrid)
      if(is.null(user.gridsize)) user.gridsize = tail(observationgrid,1)
      if(sum(is.na(X)) > 0){
        stop("can not compute 'dffm' if X contains internal NAs. The 'pace' method works for data with internal NAs")
      }
      lambda = 10^{-3}
      g = seq(from = head(observationgrid,1), to = tail(observationgrid,1), by = tail(observationgrid,1)/user.gridsize)
      grid = g / tail(g,1)
      breakpoints = (observationgrid-observationgrid[1])/(observationgrid[length(observationgrid)] - observationgrid[1])
      splinebasis = create.bspline.basis(rangeval=c(0,1), nbasis= 4 + length(breakpoints) - 2, norder=4, breaks=breakpoints)
      YieldFdPar = fdPar(fdobj = splinebasis, Lfdobj = 2, lambda = lambda)
      YieldFd = smooth.basis(breakpoints, t(X), YieldFdPar)$fd
      YieldPCA = pca.fd(YieldFd, nharm = 4 + length(breakpoints) - 2, harmfdPar = fdPar(YieldFd), centerfns = TRUE)
      Fhat = YieldPCA$scores
      eigenvalues = YieldPCA$values[1:length(observationgrid)]
    }
    # pace
    if(class(X)[1] == "xts" & pace == TRUE | class(X)[1] == "mts" & pace == TRUE){
      # preliminary
      if(is.null(observationgrid))  observationgrid = as.numeric(colnames(X))
      if(length(observationgrid) != ncol(X) ) stop('please specify a valid observationgrid')
      if(is.na(as.numeric(observationgrid))[1])stop("column names have to be numbers or you have to specify a observationgrid")
      observationgrid = as.numeric(observationgrid)
      if(is.null(user.gridsize)) user.gridsize = tail(observationgrid,1)-head(observationgrid,1)+1
      K = length(observationgrid)
      Ly = lapply(seq(length.out = nrow(X)), function(i) na.omit(X[i,]))
      Lt = lapply(seq(length.out = nrow(X)), function(i) observationgrid[!is.na(X[i,])])
      if(user.gridsize< K + 2) user.gridsize = K+2
      if(is.null(user.gridsize)){
        f = FPCA(Ly, Lt, list(error=FALSE, usergrid=TRUE, methodXi = "CE", methodSelectK=K))
      } else {
        f = FPCA(Ly, Lt, list(error=FALSE, usergrid=FALSE, methodXi = "CE", methodSelectK=K, nRegGrid=user.gridsize))
      }
      Fhat = f$xiEst
      eigenvalues = f$lambda
    }
    ## data.frame
    # fd
    if(class(X)[1] == "data.frame" & pace == FALSE){
      data = ts(X)
      data = ts_xts(data)
      if(is.null(observationgrid))  observationgrid = as.numeric(colnames(data))
      if(length(observationgrid) != ncol(data)) stop('please specify a valid observationgrid')
      if(is.na(as.numeric(observationgrid))[1])stop("column names have to be numbers or you have to specify a observationgrid")
      observationgrid = as.numeric(observationgrid)
      if(is.null(user.gridsize)) user.gridsize = tail(observationgrid,1)
      if(sum(is.na(X)) > 0) data = na.omit(data)
      lambda = 10^{-3}
      g = seq(from = head(observationgrid,1), to = tail(observationgrid,1), by = tail(observationgrid,1)/user.gridsize)
      grid = g / tail(g,1)
      breakpoints = (observationgrid-observationgrid[1])/(observationgrid[length(observationgrid)] - observationgrid[1])
      splinebasis = create.bspline.basis(rangeval=c(0,1), nbasis= 4 + length(breakpoints) - 2, norder=4, breaks=breakpoints)
      YieldFdPar = fdPar(fdobj = splinebasis, Lfdobj = 2, lambda = lambda)
      YieldFd = smooth.basis(breakpoints, t(data), YieldFdPar)$fd
      YieldPCA = pca.fd(YieldFd, nharm = 4 + length(breakpoints) - 2, harmfdPar = fdPar(YieldFd), centerfns = TRUE)
      Fhat = YieldPCA$scores
      eigenvalues = YieldPCA$values[1:length(observationgrid)]
    }
    # pace
    if(class(X)[1] == "data.frame" & pace == TRUE){
      data = ts(X)
      if(is.null(observationgrid))  observationgrid = as.numeric(colnames(data))
      if(length(observationgrid) != ncol(data) ) stop('please specify a valid observationgrid')
      if(is.na(as.numeric(observationgrid))[1])stop("column names have to be numbers or you have to specify a observationgrid")
      observationgrid = as.numeric(observationgrid)
      if(is.null(user.gridsize)) user.gridsize = tail(observationgrid,1)-head(observationgrid,1)+1
      K = length(observationgrid)
      Ly = lapply(seq(length.out = nrow(data)), function(i) na.omit(data[i,]))
      Lt = lapply(seq(length.out = nrow(data)), function(i) observationgrid[!is.na(data[i,])])
      if(is.null(user.gridsize)){
        f = FPCA(Ly, Lt, list(error=FALSE, usergrid=TRUE, methodXi = "CE", methodSelectK=K))
      } else {
        f = FPCA(Ly, Lt, list(error=FALSE, usergrid=FALSE, methodXi = "CE", methodSelectK=K, nRegGrid=user.gridsize))
      }
      Fhat = f$xiEst
      eigenvalues = f$lambda
    }
    ## Evaluate IC
    T = dim(Fhat)[1]
    q = dim(Fhat)[1]
    if(q<100)q=100
    K.max = dim(Fhat)[2]
    p.max= 8*floor(q/100)^(1/4) # according to Schwert (1989)
    eval.criterion = function(L,m){
      if(L*m+10 > T) (return(Inf))
      factors = Fhat[,1:L, drop=FALSE]
      colnames(factors) = paste("F.",1:L, sep = "")
      if(L == 1){
        aux = ar(factors, aic=FALSE, order.max = m, method = "ols", demean=FALSE)
        tr.sigmahat = sum((aux$resid)^2,na.rm=TRUE)/(T-(L*m))
      } else {
        aux=vars::VAR(factors, p=m, type="none")
        tr.sigmahat = sum((residuals(aux))^2)/(T-(L*m))
      }
      tot.prederror = log(tr.sigmahat + sum(eigenvalues[-c(1:L)]))
      return(tot.prederror + (L*m)*log(T)/T)
    }
    IC = matrix(nrow = p.max, ncol=K.max, dimnames = list(paste("p.",1:p.max,sep=""), paste("K.",1:K.max,sep="")))
    for(i in 1:p.max) (IC[i,] = sapply(1:K.max, eval.criterion, m=i))
    IC.min = which(IC == min(na.omit(IC)), arr.ind = TRUE)[1,]
    p = IC.min[1]
    K = as.numeric(IC.min[2])
  }
  ## fd object
  if(class(X)[1] == "fd"){
    # preliminary
    if(pace == TRUE & class(X)[1] == "fd")warning("input object already had been an object of class 'fd'")
    if(!is.null(observationgrid))stop("please change a new observationgrid in dffm.preprocessing")
    if(!is.null(user.gridsize))stop("please change a new user grid size in dffm.preprocessing")
    observationgrid = as.numeric(X$fdnames$time)
    user.gridsize = tail(observationgrid,1)
    if(is.null(K)) K = 3
    if(sum(is.na(X$coefs))>0){
      stop("can not compute dffm because of NA in the 'fd' object. Please use the pace method")}
    data = NULL
    # grid
    g = seq(from = head(observationgrid,1), to = tail(observationgrid,1), by = tail(observationgrid,1)/user.gridsize)
    grid = g / tail(g,1)
    YieldPCA = pca.fd(X, nharm = X$basis$nbasis, harmfdPar = fdPar(X), centerfns = TRUE)
    # output objects
    eigenvalues = YieldPCA$values[1:K]
    eigenfunctions <- eval.fd(grid, YieldPCA$harmonics)[,1:K]
    meanfunction <- eval.fd(grid, YieldPCA$meanfd)
    factorscores <- YieldPCA$scores[,1:K]
    temp <- factorscores %*% t(eigenfunctions)
    fitted.values <- matrix(NA, ncol = dim(meanfunction)[1], nrow = dim(X$coefs)[2])
    for (i in 1:dim(X$coefs)[2]){
      fitted.values[i,] <- temp[i,] + meanfunction}
    # making output objects interpretable
    if(K == 1){
      factorscores = as.matrix(factorscores)
      eigenfunctions = as.matrix(eigenfunctions)
      colnames(eigenfunctions) = "PC1"
    }
    row.names(eigenfunctions) = g
    row.names(meanfunction) = g
    dimnames(fitted.values) = list(1:dim(X$coefs)[2], g)
    dimnames(factorscores) = list(1:dim(X$coefs)[2], paste("Factor", 1:K))
  }
  ## FPCA object
  if(class(X)[1] == "FPCA"){
    # preliminary
    if(!is.null(observationgrid))stop("please change a new observationgrid in dffm.preprocessing")
    if(!is.null(user.gridsize))stop("please change a new user grid size in dffm.preprocessing")
    observationgrid = X$obsGrid
    if(length(observationgrid) != length(X$obsGrid)) stop('please specify a valid observationgrid')
    if(is.null(K)) K = 3
    data = NULL
    g = X$workGrid
    # output objects
    eigenvalues = X$lambda[1:K]
    eigenfunctions = X$phi[,1:K]
    meanfunction = X$mu
    factorscores = X$xiEst[,1:K]
    fitted.values = matrix(data = rep(meanfunction, nrow(X$xiEst)), nrow = nrow(X$xiEst),
                           byrow = TRUE) + factorscores %*% t(eigenfunctions)
    # making output objects interpretable
    if(K == 1){
      eigenfunctions = as.matrix(eigenfunctions)
      factorscores = as.matrix(factorscores)
    }
    dimnames(eigenfunctions) = list(g, paste("PC", 1:K))
    meanfunction = matrix(meanfunction, dimnames = list(g, "mean"))
    dimnames(factorscores) = list(1:dim(X$xiEst)[1], paste("Factor", 1:K))
    dimnames(fitted.values) = list(1:dim(X$xiEst)[1], g)
  }
  ## ts/xts object
  # fd
  if(class(X)[1] == "xts" & pace == FALSE | class(X)[1] == "mts" & pace == FALSE){
    # preliminary
    data = X
    if(is.null(observationgrid))  observationgrid = as.numeric(colnames(data))
    if(length(observationgrid) != ncol(data) ) stop('please specify a valid observationgrid')
    if(is.na(as.numeric(observationgrid))[1])stop("column names have to be numbers or you have to specify a observationgrid")
    observationgrid = as.numeric(observationgrid)
    if(is.null(user.gridsize)) user.gridsize = tail(observationgrid,1)
    if(sum(is.na(X)) > 0){
      stop("can not compute 'dffm' if X contains internal NAs. The 'pace' method works for data with internal NAs")
    }
    lambda = 10^{-3}
    if(is.null(K)) K = 3
    g = seq(from = head(observationgrid,1), to = tail(observationgrid,1), by = tail(observationgrid,1)/user.gridsize)
    grid = g / tail(g,1)
    breakpoints = (observationgrid-observationgrid[1])/(observationgrid[length(observationgrid)] - observationgrid[1])
    splinebasis = create.bspline.basis(rangeval=c(0,1), nbasis= 4 + length(breakpoints) - 2, norder=4, breaks=breakpoints)
    YieldFdPar = fdPar(fdobj = splinebasis, Lfdobj = 2, lambda = lambda)
    YieldFd = smooth.basis(breakpoints, t(data), YieldFdPar)$fd
    YieldPCA = pca.fd(YieldFd, nharm = 4 + length(breakpoints) - 2, harmfdPar = fdPar(YieldFd), centerfns = TRUE)
    # output objects
    eigenvalues = YieldPCA$values[1:K]
    eigenfunctions <- eval.fd(grid, YieldPCA$harmonics)[,1:K]
    meanfunction <- eval.fd(grid, YieldPCA$meanfd)
    factorscores <- YieldPCA$scores[,1:K]
    temp <- factorscores %*% t(eigenfunctions)
    fitted_yields <- matrix(NA, ncol = dim(meanfunction)[1], nrow = dim(data)[1])
    for (i in 1:dim(data)[1]){
      fitted_yields[i,] <- temp[i,] + meanfunction
    }
    # making output objects interpretable
    if(K == 1){
      factorscores = as.matrix(factorscores)
      eigenfunctions = as.matrix(eigenfunctions)
      colnames(eigenfunctions) = "PC1"
    }
    colnames(factorscores) = paste("Factor", 1:K)
    factorscores = ts(factorscores, start = head(time(data),1), frequency = frequency(data))
    data = ts_xts(X)
    output.y = xts(x = fitted_yields, order.by = time(data))
    row.names(eigenfunctions) = g
    row.names(meanfunction) = g
    colnames(output.y) = g
    fitted.values = ts_ts(output.y)
    data = ts_ts(data)
  }
  # pace
  if(class(X)[1] == "xts" & pace == TRUE | class(X)[1] == "mts" & pace == TRUE){
    # preliminary
    data = ts_ts(X)
    if(is.null(observationgrid))  observationgrid = as.numeric(colnames(data))
    if(length(observationgrid) != ncol(data) ) stop('please specify a valid observationgrid')
    if(is.na(as.numeric(observationgrid))[1])stop("column names have to be numbers or you have to specify a observationgrid")
    observationgrid = as.numeric(observationgrid)
    if(is.null(user.gridsize)) user.gridsize = tail(observationgrid,1)-head(observationgrid,1)+1
    if(is.null(K)) K = 3
    if(user.gridsize< K + 2) user.gridsize = K+2
    Ly = lapply(seq(length.out = nrow(X)), function(i) na.omit(X[i,]))
    Lt = lapply(seq(length.out = nrow(X)), function(i) observationgrid[!is.na(X[i,])])
    if(is.null(user.gridsize)){
      f = FPCA(Ly, Lt, list(error=FALSE, usergrid=TRUE, methodXi = "CE", methodSelectK=K))
    } else {
      f = FPCA(Ly, Lt, list(error=FALSE, usergrid=FALSE, methodXi = "CE", methodSelectK=K, nRegGrid=user.gridsize))
    }
    # output objects
    workGrid = round(f$workGrid,2)
    g = workGrid
    eigenvalues = f$lambda
    eigenfunctions = f$phi
    meanfunction = matrix(f$mu, dimnames = list(workGrid,"mean"))
    factorscores = f$xiEst
    fitted.values = matrix(rep(meanfunction, nrow(data)), nrow=nrow(data), byrow=TRUE) + factorscores %*% t(eigenfunctions)
    # making output objects interpretable
    dimnames(eigenfunctions) = list(workGrid, paste("PC",1:K,sep=""))
    colnames(factorscores) = paste("Factor", 1:K)
    factorscores = ts(factorscores, start = head(time(data),1), frequency = frequency(data))
    data = ts_xts(X)
    output.y = xts(x = fitted.values, order.by = time(data))
    colnames(output.y) = g
    fitted.values = ts_ts(output.y)
    data = ts_ts(data)
  }
  ## data.frame object
  # fd
  if(class(X)[1] == "data.frame" & pace == FALSE){
    # preliminary
    data = ts(X)
    data = ts_xts(data)
    if(is.null(observationgrid))  observationgrid = as.numeric(colnames(data))
    if(length(observationgrid) != ncol(data) ) stop('please specify a valid observationgrid')
    if(is.na(as.numeric(observationgrid))[1])stop("column names have to be numbers or you have to specify a observationgrid")
    observationgrid = as.numeric(observationgrid)
    if(is.null(user.gridsize)) user.gridsize = tail(observationgrid,1)
    if(sum(is.na(X)) > 0){
      data = na.omit(data)
      warning("missing data in X had been removed")
    }
    lambda = 10^{-3}
    if(is.null(K)) K = 3
    g = seq(from = head(observationgrid,1), to = tail(observationgrid,1), by = tail(observationgrid,1)/user.gridsize)
    grid = g / tail(g,1)
    breakpoints = (observationgrid-observationgrid[1])/(observationgrid[length(observationgrid)] - observationgrid[1])
    splinebasis = create.bspline.basis(rangeval=c(0,1), nbasis= 4 + length(breakpoints) - 2, norder=4, breaks=breakpoints)
    YieldFdPar = fdPar(fdobj = splinebasis, Lfdobj = 2, lambda = lambda)
    YieldFd = smooth.basis(breakpoints, t(data), YieldFdPar)$fd
    YieldPCA = pca.fd(YieldFd, nharm = 4 + length(breakpoints) - 2, harmfdPar = fdPar(YieldFd), centerfns = TRUE)
    # output objects
    eigenvalues = YieldPCA$values[1:K]
    eigenfunctions <- eval.fd(grid, YieldPCA$harmonics)[,1:K]
    meanfunction <- eval.fd(grid, YieldPCA$meanfd)
    factorscores <- YieldPCA$scores[,1:K]
    temp <- factorscores %*% t(eigenfunctions)
    fitted.values <- matrix(NA, ncol = dim(meanfunction)[1], nrow = dim(data)[1])
    for (i in 1:dim(data)[1]){
      fitted.values[i,] <- temp[i,] + meanfunction
    }
    # making output objects interpretable
    if(K == 1){
      factorscores = as.matrix(factorscores)
      eigenfunctions = as.matrix(eigenfunctions)
    }
    dimnames(eigenfunctions) = list(g, paste("PC",1:K,sep=""))
    row.names(meanfunction) = g
    dimnames(factorscores) = list(1:dim(data)[1], paste("Factor", 1:K))
    dimnames(fitted.values) = list(1:dim(data)[1], g)
    data = as.data.frame(data)
  }
  # pace
  if(class(X)[1] == "data.frame" & pace == TRUE){
    # preliminary
    data = ts(X)
    if(is.null(observationgrid))  observationgrid = as.numeric(colnames(data))
    if(length(observationgrid) != ncol(data) ) stop('please specify a valid observationgrid')
    if(is.na(as.numeric(observationgrid))[1])stop("column names have to be numbers or you have to specify a observationgrid")
    observationgrid = as.numeric(observationgrid)
    if(is.null(user.gridsize)) user.gridsize = tail(observationgrid,1)-head(observationgrid,1)+1
    if(is.null(K)) K = 3
    Ly = lapply(seq(length.out = nrow(data)), function(i) na.omit(data[i,]))
    Lt = lapply(seq(length.out = nrow(data)), function(i) observationgrid[!is.na(data[i,])])
    if(is.null(user.gridsize)){
      f = FPCA(Ly, Lt, list(error=FALSE, usergrid=TRUE, methodXi = "CE", methodSelectK=K))
    } else {
      f = FPCA(Ly, Lt, list(error=FALSE, usergrid=FALSE, methodXi = "CE", methodSelectK=K, nRegGrid=user.gridsize))
    }
    # output objects
    workGrid = round(f$workGrid,2)
    g = workGrid
    eigenvalues = f$lambda
    eigenfunctions = f$phi
    meanfunction = matrix(f$mu, dimnames = list(workGrid,"mean"))
    factorscores = f$xiEst
    fitted.values = matrix(rep(meanfunction, nrow(data)), nrow=nrow(data), byrow=TRUE) + factorscores %*% t(eigenfunctions)
    # making output objects interpretable
    dimnames(eigenfunctions) = list(g, paste("PC", 1:K))
    row.names(meanfunction) = g
    dimnames(factorscores) = list(1:dim(data)[1], paste("Factor", 1:K))
    dimnames(fitted.values) = list(1:dim(data)[1], g)
    data = as.data.frame(data)
  }
  ## output
  if(is.null(obse.name)) obse.name = "NULL"
  if(is.null(grid.name)) grid.name = "NULL"
  if(is.null(K.name)) K.name = "NULL"
  if(!is.null(p)) p = as.numeric(p)
  string = "dffm(X = %s,pace = %s,criterion = %s,observationgrid = %s,user.gridsize = %s,K = %s)"
  vals = c(as.character(X.name),as.character(pace), as.character(criterion.name),as.character(obse.name), as.character(grid.name), as.character(K.name))
  call = noquote(do.call(sprintf, as.list(c(string, vals))))
  basis = list("call" = call, "observationgrid" = observationgrid, "workgrid" = g, "p.opt" = p, "K" =  K, "data" = data)
  # output list
  output = list(
    "basis" = basis,
    "eigenvalues" = eigenvalues,
    "loadingfunctions" = eigenfunctions,
    "meanfunction" = meanfunction,
    "factorscores" = factorscores,
    "fitted.data" = fitted.values)
  class(output) = "dffm"
  return(invisible(output))
}

#' print dffm
#'
#' The generic S3 method print for an object of class 'dffm'.
#'
#' @param x An object which has to be of class 'dffm'.
#'
#' @return
#' A comprehension, which shows the class and dimensions of the list objects, stored in the 'dffm' object.
#' @export
#' @examples
#' data = load.dieboldli()
#' d = dffm(data)
#' print(d)
#' # or
#' d
print.dffm = function(x){
  print(list(
    "eigenvalues" = paste("class:", class(x$eigenvalues)[1], ";", "length:",length(x$eigenvalues)),
    "loadingfunctions" = paste("class:", class(x$loadingfunctions)[1],";",dim(x$loadingfunctions)[1],
                             "rows,", dim(x$loadingfunctions)[2], "cols"),
    "meanfunction" = paste("class:", class(x$meanfunction)[1],";",dim(x$meanfunction)[1],
                           "rows,", dim(x$meanfunction)[2], "cols"),
    "factorscores" = paste("class:", class(x$factorscores)[1],";",dim(x$factorscores)[1], "rows,",
                           dim(x$factorscores)[2], "cols"),
    "fitted.data" = paste("class:", class(x$fitted.data)[1],";",dim(x$fitted.data)[1], "rows,",
                          dim(x$fitted.data)[2], "cols")
  ), quote = FALSE)
}

#' summary dffm
#'
#' The generic S3 method summary for an object of class 'dffm'.
#'
#' @param x An object which has to be of class 'dffm'.
#'
#' @return
#' Will print 'eigenvalues', 'call', 'observationgrid', 'workinggrid', 'K', and 'p.opt' of the 'dffm' object.
#' @export
#' @examples
#' data = load.dieboldli()
#' d = dffm(data)
#' summary(d)
summary.dffm = function(x){
  output = list("eigenvalues" = x$eigenvalues, "call" = x$basis$call, "observationgrid" = x$basis$observationgrid,
                "workinggrid" = x$basis$workgrid, "K" = x$basis$K, "p.opt" = x$basis$p.opt)
  return(output)
}

#' plot dffm
#'
#' The generic S3 method plot for an object of class 'dffm'.
#'
#' @param x An object which has to be of class 'dffm'.
#' @param time A parameter which selects the wanted fitted data observation. If time is NULL, time will be 1. NULL is default
#' @param date.on Will blend in the date of the observation, if it was from a time series object. If date.on is FALSE
#' the plot will not show the date. FALSE is default.
#' @param ... Further plot inputs, such as 'ylab' or 'main' are possible.
#'
#' @return
#' A plot of the fitted data.
#' @export
#' @examples
#' data = load.dieboldli()
#' d = dffm(data)
#' plot(d, time = 20, date.on = TRUE)
#' # or
#' plot(d, time = 5, xlab = "maturities", ylab = "yields", main = "yield curve of mai 1970")
plot.dffm = function(x, time = NULL, date.on = FALSE, ...){
  if(is.null(time)) time = 1
  if(date.on == TRUE & class(x$fitted.data)[1] != "mts") stop("'dffm' object has no time series 'fitted.data'")
  grid = as.numeric(x$basis$workgrid)
  y = x$fitted.data[time,]
  plot(x = grid, y = y, type = "l")

  if(class(x$fitted.data)[1] == "mts" & date.on == TRUE){
    plot(x = grid, y = y, type = "l",
         sub = list(paste("from", as.yearmon(time(x$fitted.data))[time]), cex =.9), ...)
  }
  if(class(x$fitted.data)[1] == "mts" & date.on == FALSE){
    plot(x = grid, y = y, type = "l", ...)
  }
  if(class(x$fitted.data)[1] != "mts"){
    plot(x = grid, y = y, type = "l", ...)
  }
}

#' dffm preprocessing
#'
#' dffm preprocessing is used to prepare a time series or data frame data object especially for the dffm function. For more details
#' see dffm. The output objects are based on the fda package ('fd') or on the fdapace package ('FPCA').
#'
#' @param data A multivariate time series of class 'ts' or 'data.frame'.
#' @param pace Preprocessing method to be used. If TRUE the PACE method from the package 'fdapace' is applied (see Yao et al., 2005).
#' In case of irregularly observed functional data please set to TRUE.
#' If FALSE smoothing splines from the package 'fda' are used (see Kokoszka and Reimherr, 2017). FALSE is default.
#' @param observationgrid An optional list, sequence or numeric vector which contains a predetermined observationgrid for the
#' observations. If NULL the column names will be used. NULL is default.
#' @param user.gridsize An optional parameter specifying the size of the workinggrid. If user.gridsize is NULL, user.gridsize will be the highest
#' number from the column names or the specified observationgrid. NULL is default.
#' @param lambda A smoothing parameter for the penalty term which prevents over fitting. If lambda is 0, the penalty term will
#' not contribute to the analysis. Lambda is only used in the 'FUNCTIONAL DATA ANALYSIS' ('fd') method. If lambda is NULL,
#' lambda will be 10^{-3}. NULL is default.
#'
#' @return
#' An object of class 'fd' if pace is FALSE, otherwise an object of class 'FPCA'.
#' @export
#' @references
#' * Kokoszka P. and Reimherr M. (2017), "Introduction to Functional Data Analysis". CRC Press.
#' * Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). Functional data analysis for sparse longitudinal
#' data. Journal of the American Statistical Association, 100:577-590.
#' @examples
#' # standard case
#' dl = load.dieboldli()
#' pp = dffm.preprocessing(data = dl, pace = FALSE, observationgrid = NULL, lambda = 10^{-4})
#' x = dffm(pp, K = 3)
#' # pace method
#' pp = dffm.preprocessing(data = dl, pace = TRUE, user.gridsize = 180)
#' x = dffm(pp)
dffm.preprocessing = function(data, pace = FALSE, observationgrid = NULL, user.gridsize = NULL, lambda = NULL){
  if(missing(data)) stop("please provide valid data")
  if(pace==FALSE)if(!is.null(observationgrid))if(length(observationgrid)==ncol(data))colnames(data)=as.numeric(observationgrid)
  if(is.null(observationgrid))  observationgrid = as.numeric(colnames(data))
  if(length(observationgrid) != ncol(data)) stop('please specify a valid observationgrid')
  if(is.na(as.numeric(observationgrid))[1])stop("column names have to be numbers or you have to specify a observationgrid")
  observationgrid = as.numeric(observationgrid)
  if(pace == FALSE) if(is.null(user.gridsize)) user.gridsize = tail(observationgrid,1)
  if(pace == TRUE) if(is.null(user.gridsize)) user.gridsize = tail(observationgrid,1)-head(observationgrid,1)+1
  if(!is.null(lambda) & pace == TRUE)warning("there is no lambda to specify for the 'FPCA' object")
  if(pace == FALSE & sum(is.na(data))>0){
    warning("because of missing data in 'data', the 'fd' object will have NAs. The 'pace' method will not omit NAs")}
  # fd
  if(is.null(lambda)) lambda = 10^{-3}
  breakpoints = (observationgrid-observationgrid[1])/(observationgrid[length(observationgrid)] - observationgrid[1])
  splinebasis = create.bspline.basis(rangeval=c(0,1), nbasis= 4 + length(breakpoints) - 2, norder=4, breaks=breakpoints)
  YieldFdPar = fdPar(fdobj = splinebasis, Lfdobj = 2, lambda = lambda)
  output = smooth.basis(breakpoints, t(data), YieldFdPar)$fd
  # pace
  K = length(observationgrid)
  if(pace == TRUE){
    if(class(data)[1] == "data.frame") data = ts(data)
    Ly = lapply(seq(length.out = nrow(data)), function(i) na.omit(data[i,]))
    Lt = lapply(seq(length.out = nrow(data)), function(i) observationgrid[!is.na(data[i,])])
    if(is.null(user.gridsize)){
      f = FPCA(Ly, Lt, list(error=FALSE, usergrid=TRUE, methodXi = "CE", methodSelectK=K))
    }else{
      f = FPCA(Ly, Lt, list(error=FALSE, usergrid=FALSE, methodXi = "CE", methodSelectK=K, nRegGrid=user.gridsize))
    }
    output = f
  }
  return(output)
}

#' dffm criterion
#'
#' Computes an optimal number of factors K and an optimal number of lags p for functional data analysis or prediction, based on the
#' information criterion by Otto and Salish (2021).
#'
#' @param x An object which has to be of class 'dffm'.
#' @param K.max A predetermined maximum of K for the criterion.
#' @param p.max A predetermined maximum of p for the criterion.
#'
#' @return
#' An optimal number of K and p based on the information criteria from Otto and Salish (2021) and all possible combinations of K
#' and p stored in a matrix.
#'
#' @export
#' @references
#' * Otto S. and Salish N. (2021). "Dynamic Factor Model for Functional Time Series: Identification, Estimation, and Prediction". ...
#' @examples
#' data = load.dieboldli()
#' d = dffm(d)
#' dffm.criterion(d)
dffm.criterion = function(x, K.max=NULL, p.max=NULL){
  if(class(x) != "dffm") stop('x must be class "dffm"')
  Fhat = x$factorscores
  T = dim(Fhat)[1]
  q = dim(Fhat)[1]
  if(q<100)q=100
  if(is.null(K.max)){
    K.max = dim(Fhat)[2]
  }
  if(is.null(p.max))  p.max= 8*floor(q/100)^(1/4) # according to Schwert (1989)
  ## Evaluate IC
  eval.criterion = function(L,m){
    if(L*m+10 > T) (return(Inf))
    factors = Fhat[,1:L, drop=FALSE]
    colnames(factors) = paste("F.",1:L, sep = "")
    if(L == 1){
      aux = ar(factors, aic=FALSE, order.max = m, method = "ols", demean=FALSE)
      tr.sigmahat = sum((aux$resid)^2,na.rm=TRUE)/(T-(L*m))
    } else {
      aux=vars::VAR(factors, p=m, type="none")
      tr.sigmahat = sum((residuals(aux))^2)/(T-(L*m))
    }
    tot.prederror = log(tr.sigmahat + sum(x$eigenvalues[-c(1:L)]))
    return(tot.prederror + (L*m)*log(T)/T) # (i) L*m (ii) L+m, (iii) L+m+L*m
  }
  IC = matrix(nrow = p.max, ncol=K.max, dimnames = list(paste("p.",1:p.max,sep=""), paste("K.",1:K.max,sep="")))
  for(i in 1:p.max) (IC[i,] = sapply(1:K.max, eval.criterion, m=i))
  IC.min = which(IC == min(na.omit(IC)), arr.ind = TRUE)[1,]
  names(IC.min) = c("p.opt", "K.opt")
  return(list("IC.min" = IC.min,"IC" = IC))
}

VAR.forecast = function(x, AR = FALSE, criterion = FALSE, h = NULL, p = NULL, K = NULL){
  if(class(x) != "dffm") stop('x must be class "dffm"')
  if(criterion == TRUE & !is.null(p) | criterion == TRUE & !is.null(K)| criterion == TRUE & !is.null(K) & !is.null(p)){
    stop("please use either criterion or predetermined K and/or p")}
  if(criterion == TRUE){
    IC = dffm.criterion(x)
    p = as.numeric(IC$IC.min[1])
    K = as.numeric(IC$IC.min[2])
  }
  if(criterion == FALSE & is.null(K)){K = length(x$eigenvalues)}
  if(criterion == FALSE & is.null(p)){p = 8}
  if(is.null(h)){h = 1}
  if(K > length(x$eigenvalues)) stop("K is greater than number of components in the 'dffm' object")
  factors = x$factorscores[,1:K]
  if(K == 1) factors = as.matrix(factors)
  if(AR==FALSE){
    var.model<-ar(factors, aic = FALSE, order.max = p, method = "ols")
    predicted.factors = suppressWarnings(predict(var.model, factors, n.ahead = h)$pred)
  } else {
    sep.pred = list()
    for(l in 1:K){
      ar.model<-ar(factors[,l], aic = FALSE, order.max = p, method = "ols")
      sep.pred[[l]] = predict(ar.model, factors[,l], n.ahead = h)$pred
    }
    predicted.factors = do.call(cbind, sep.pred)
  }
  if(class(x$fitted.data)[1] == "mts"){
    predicted.factors = ts(predicted.factors, start = tail(time(x$fitted.data),1) + 1/frequency(x$fitted.data), frequency = frequency(x$fitted.data))
  }
  if(class(x$fitted.data)[1] != "mts"){rownames(predicted.factors) = (dim(x$fitted.data)[1]+1):(dim(x$fitted.data)[1]+dim(temp)[1])}
  return(predicted.factors)
}

#' dffm predict
#'
#' Will make predictions of the fitted input data from a 'dffm' object based on vector autoregressive or autoregressive models.
#'
#' @param x An object which has to be of class 'dffm'.
#' @param AR Will determine if prediction is based on vector autoregressive models or autoregressive models. If AR is FALSE
#' prediction will be based on vector autoregressive models. FALSE is default.
#' @param criterion Automatic selection of the number of factors K and lags p. If TRUE the consistent information criterion of Otto and Salish (2021)
#' is used. If FALSE the value specified in 'K' and 'p' is used. FALSE is default. For more details see 'dffm.criterion'.
#' @param h Selects the amount of predictions which will be made. If h is NULL, h will be 1. NULL is default.
#' @param p Selects the lags which will be used in the prediction. If p is NULL, p will be 8. NULL is default. If criterion is
#' TRUE p will be overwritten.
#' @param K Selects the number of factors which will be used for the prediction. If K is NULL, K will be maximum of possible
#' components. NULL is default. If criterion is TRUE K will be overwritten.
#'
#' @return
#' \item{predicted.values}{A matrix or a ts-matrix containing the predicted fitted data.}
#' \item{predicted.factors}{A matrix or a ts-matrix containing the predicted factors.}
#' @export
#' @references
#' * Otto S. and Salish N. (2021). "Dynamic Factor Model for Functional Time Series: Identification, Estimation, and Prediction". ...
#' @examples
#' data = load.dieboldli()
#' d = dffm(d)
#' dffm.predict(d)
dffm.predict = function(x, AR = FALSE, criterion = FALSE, h = NULL, p = NULL, K = NULL){
  if(class(x) != "dffm") stop('x must be class "dffm"')
  if(criterion == TRUE & !is.null(p) | criterion == TRUE & !is.null(K)| criterion == TRUE & !is.null(K) & !is.null(p)){
    stop("please use either criterion or predetermined K and/or p")}
  if(criterion == TRUE){
    IC = dffm.criterion(x)
    p = as.numeric(IC$IC.min[1])
    K = as.numeric(IC$IC.min[2])
  }
  if(K > length(x$eigenvalues)) stop("K is greater then the number of possible components")
  if(K < 0) stop("K is not allowed to be negative")
  if(class(x$fitted.data)[1] == "mts"){
    if(sum(is.na(x$fitted.data))>0){
      x$factorscores = ts_xts(x$factorscores)
      x$factorscores = na.omit(x$factorscores)
      x$fitted.data = ts_xts(x$fitted.data)
      x$fitted.data = na.omit(x$fitted.data)}
  }
  if(class(x$fitted.data)[1] != "mts"){
    x$factorscores = na.omit(x$factorscores)
    x$fitted.data = na.omit(x$fitted.data)
  }
  if(criterion == FALSE & is.null(K)){K = length(x$eigenvalues)}
  if(criterion == FALSE & is.null(p)){p = 8}
  if(is.null(h)){h = 1}
  f.cast = VAR.forecast(x, h = h, AR = AR, p = p, K = K)
  temp = f.cast %*% t(x$loadingfunctions[,1:K])
  fitt = matrix(NA, ncol = dim(temp)[2], nrow = dim(temp)[1])
  for (i in 1:dim(temp)[1]) {
    fitt[i,] = temp[i,] + x$meanfunction
  }
  if(class(x$fitted.data)[1] == "mts"){
    fitt = ts(fitt, start = tail(time(x$fitted.data),1) + 1/frequency(x$fitted.data), frequency = frequency(x$fitted.data))
  }
  if(class(x$fitted.data)[1] != "mts"){rownames(fitt) = (dim(x$fitted.data)[1]+1):(dim(x$fitted.data)[1]+dim(temp)[1])}
  if(K == 0){
    fitt = matrix(NA, nrow = h, ncol = dim(x$meanfunction)[1])
    for (i in 1:h){
      fitt[i,] = x$meanfunction
    }
    if(class(x$fitted.data)[1] == "mts"){
      fitt = ts(fitt, start = tail(time(x$fitted.data),1) + 1/frequency(x$fitted.data), frequency = frequency(x$fitted.data))}
    if(class(x$fitted.data)[1] != "mts"){rownames(fitt) = (dim(x$fitted.data)[1]+1):(dim(x$fitted.data)[1]+dim(fitt)[1])}}
  colnames(fitt) = x$basis$workgrid
  fitted.values = fitt
  output = list("predicted.values" = fitted.values, "predicted.factors" = f.cast)
  return(output)
}


#' load dieboldli
#'
#' A function which takes yield curve data provided by Diebold and Li (2006): https://www.sas.upenn.edu/~fdiebold/papers/paper49/FBFITTED.txt
#'
#' @return
#' A time series of yield curve data with maturities ranging from 1 to 120 month.
#' @references
#' * Diebold, F. X. and Li, C. (2006). Forecasting the term structure of government bond yields. Journal of Econometrics, 130:337-364.
#' @export
#' @examples
#' data = load.dieboldli()
load.dieboldli = function(){
  data <- read.csv("https://www.sas.upenn.edu/~fdiebold/papers/paper49/FBFITTED.txt", skip =12, header=FALSE, sep="")
  yields<- matrix(na.omit(as.vector(t(cbind(data[c(TRUE,FALSE),2:13], data[c(FALSE,TRUE),1:7])))), ncol=18, byrow=TRUE)
  Y <- ts(yields, start=c(1970, 1), frequency=12)
  colnames(Y) = c(1,3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
  return(Y)
}

#' load jungbacker
#'
#' A function which takes yield curve data provided by Jungbacker et.al (2014): http://qed.econ.queensu.ca/jae/2014-v29.1/jungbacker-koopman-van_der_wel/
#'
#' @return
#' A time series of yield curve data with maturities ranging from 3 to 120 month.
#' @references
#' * Jungbacker, B., Koopman, S. J., and Van der Wel, M. (2014). Smooth dynamic factor analysis with application to the US term structure of
#' interest rates. Journal of Applied Econometrics, 29:65-90.
#' @export
#' @examples
#' data = load.jungbacker()
load.jungbacker = function(){
  temp <- tempfile()
  download.file("http://qed.econ.queensu.ca/jae/2014-v29.1/jungbacker-koopman-van_der_wel/jkv-data.zip",temp)
  data <- read.table(unz(temp, "UnsmFB_70-09.txt"))
  unlink(temp)
  jungbacker.data = ts(data, start = c(1970, 1), frequency = 12)
  colnames(jungbacker.data) = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
  return(jungbacker.data)
}
