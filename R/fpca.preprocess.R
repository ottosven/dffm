#' fpca.preprocess
#'
#' fpca.preprocess is used to prepare a time series or data frame data object especially for the dffm function. For more details
#' see dffm. The output objects are based on the fda package or on the fdapace package.
#'
#' @param data A multivariate time series of class 'ts' or 'data.frame'.
#' @param workinggrid An optional list, sequence or numeric vector which contains a predetermined workinggrid for the analysis.
#' If NULL a sequence, beginning from the smallest to the highest observed observation will be used. NULL is default.
#' @param observationgrid An optional list, sequence or numeric vector which contains a predetermined observationgrid for the
#' observations. If NULL the column names will be used. NULL is default.
#' @param method Preprocessing method to be used. If 'splines' smoothing splines from the package 'fda' are used (see Kokoszka
#' and Reimherr, 2017). If 'pace' the PACE method from the package 'fdapace' is applied (see Yao et al., 2005). In case of
#' irregularly observed functional data please set to 'pace'. 'splines' is default.
#'
#' @return
#' fpca.preprocess will return an object of class 'FPCAobj' with:
#' \item{eigenvalues}{A vector containing the eigenvalues of the sample covariance operator.}
#' \item{fpcascores}{A matrix/ts-matrix of the computed factorscores.}
#' \item{eigenfunctions.obsgrid}{A matrix containing the estimated loading functions, which are evaluated on the workinggrid.}
#' \item{meanfunction.obsgrid}{A vector containing the estimated mean function, which is evalauted on the observationgrid.}
#' \item{eigenfunctions.workgrid}{A matrix containing the estimated loading functions, which are evaluated on the observationgrid.}
#' \item{meanfunction.workgrid}{A vector containing the estimated mean function, which is evalauted on the workinggrid.}
#' \item{raw.data}{Ouput of the data used in the analysis.}
#' \item{observationgrid}{Used observationgrid in the analysis.}
#' \item{workinggrid}{Used workinggrid in the analysis.}
#' \item{method}{Used method in the analysis.}
#'
#' @export
#' @references
#' * Kokoszka P. and Reimherr M. (2017), "Introduction to Functional Data Analysis". CRC Press.
#' * Otto S. and Salish N. (2022). "Dynamic Factor Model for Functional Time Series: Identification, Estimation, and Prediction". ...
#' * Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). Functional data analysis for sparse longitudinal
#' data. Journal of the American Statistical Association, 100:577-590.
#' @examples
#' # normal workinggrid
#' fed = load.fed()
#' fpca.preprocess(data = fed, method = "splines")
#' # changed workinggrid
#' fed = load.fed()
#' wg = (2:720)/2
#' fpca.preprocess(fed, workinggrid = wg, method = "splines")
fpca.preprocess = function(data, workinggrid = NULL, observationgrid = NULL, method = c("splines", "pace")){
  ## before analysis
  # checking data
  if(missing(data)) stop("please provide valid data. For more details see 'data' in the help file")
  # checking method
  if(!missing(method) & length(method) > 1) stop("Only one 'method' allowed")
  method = match.arg(method)
  # checking observationgrid
  if(is.null(observationgrid))  observationgrid = as.numeric(colnames(data))
  if(length(observationgrid) != ncol(data)) stop('please specify a valid observationgrid')
  if(is.na(as.numeric(observationgrid))[1])stop("column names have to be numbers")
  observationgrid.unit = (observationgrid-head(observationgrid,1))/(tail(observationgrid,1) - head(observationgrid,1))
  # checking workinggrid
  if(is.null(workinggrid)) workinggrid = head(observationgrid,1):tail(observationgrid,1)
  workinggrid.unit = (workinggrid-head(workinggrid,1))/(tail(workinggrid,1) - head(workinggrid,1))
  ## analysis
  # splines
  if(method == "splines"){
    splinebasis = fda::create.bspline.basis(rangeval=c(0,1), nbasis= 4 + length(observationgrid) - 2, norder=4, breaks=observationgrid.unit)
    FdPar.aux = fda::fdPar(fdobj = splinebasis)
    get.meangcv = function(lamb) (mean(fda::lambda2gcv(log10lambda = lamb, argvals = observationgrid.unit, y=t(data), fdParobj = FdPar.aux)))
    lambdas.crossval = -(1:100)/10
    meangcv = sapply(lambdas.crossval, get.meangcv)
    lambda = 10^(lambdas.crossval[which(meangcv == min(na.omit(meangcv)))])
    FdPar = fda::fdPar(fdobj = splinebasis, lambda = lambda)
    Fd = fda::smooth.basis(observationgrid.unit, t(data), FdPar)$fd
    PCA = fda::pca.fd(Fd, nharm = Fd$basis$nbasis, harmfdPar = fda::fdPar(Fd), centerfns = TRUE)
    eigenvalues = PCA$values
    fpcascores <- PCA$scores
    #
    colnames(fpcascores) = paste("PC",1:dim(fpcascores)[2], sep = "")
    if(is.ts(data)) fpcascores = ts(fpcascores, start = head(time(data),1), frequency = frequency(data))
    eigenfunctions.obsgrid <- fda::eval.fd(observationgrid.unit, PCA$harmonics)
    meanfunction.obsgrid <- fda::eval.fd(observationgrid.unit, PCA$meanfd)
    rownames(eigenfunctions.obsgrid) = observationgrid
    rownames(meanfunction.obsgrid) = observationgrid
    eigenfunctions.workgrid <- fda::eval.fd(workinggrid.unit, PCA$harmonics)
    meanfunction.workgrid <- fda::eval.fd(workinggrid.unit, PCA$meanfd)
    rownames(eigenfunctions.workgrid) = workinggrid
    rownames(meanfunction.workgrid) = workinggrid
  }
  # pace
  if(method == "pace"){
    Ly = lapply(seq(length.out = nrow(data)), function(i) na.omit(data[i,]))
    Lt = lapply(seq(length.out = nrow(data)), function(i) observationgrid[!is.na(data[i,])])
    f.obs = fdapace::FPCA(Ly, Lt, list(error=FALSE, usergrid=TRUE, methodXi = "CE", methodSelectK = dim(data)[2]+2, nRegGrid = length(workinggrid)))
    f.work = fdapace::FPCA(Ly, Lt, list(error=FALSE, usergrid=FALSE, methodXi = "CE", methodSelectK = dim(data)[2]+2, nRegGrid = length(workinggrid)))
    eigenvalues = f.work$lambda
    eigenfunctions.obsgrid = f.obs$phi
    eigenfunctions.workgrid = f.work$phi
    meanfunction.obsgrid = f.obs$mu
    meanfunction.workgrid = f.work$mu
    fpcascores = f.work$xiEst
    #
    colnames(fpcascores) = paste("PC",1:dim(fpcascores)[2], sep = "")
    if(is.ts(data)) fpcascores = ts(fpcascores, start = head(time(data),1), frequency = frequency(data))
    rownames(eigenfunctions.obsgrid) = observationgrid
    colnames(eigenfunctions.obsgrid) = paste("PC",1:dim(eigenfunctions.obsgrid)[2], sep = "")
    meanfunction.obsgrid = as.matrix(meanfunction.obsgrid)
    dimnames(meanfunction.obsgrid) = list(observationgrid, "mean")
    rownames(eigenfunctions.workgrid) = workinggrid
    colnames(eigenfunctions.workgrid) = paste("PC",1:dim(eigenfunctions.workgrid)[2], sep = "")
    meanfunction.workgrid = as.matrix(meanfunction.workgrid)
    dimnames(meanfunction.workgrid) = list(workinggrid, "mean")
  }
  ## output
  out=list(
    "eigenvalues" = eigenvalues,
    "fpcascores" = fpcascores,
    "eigenfunctions.obsgrid" = eigenfunctions.obsgrid,
    "meanfunction.obsgrid" = meanfunction.obsgrid,
    "eigenfunctions.workgrid" = eigenfunctions.workgrid,
    "meanfunction.workgrid" = meanfunction.workgrid,
    "raw.data" = data,
    "observationgrid" = observationgrid,
    "workinggrid" = workinggrid,
    "method" = method
  )
  class(out) = "FPCAobj"
  return(out)
}
