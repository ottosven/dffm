#' fpca.preprocess
#'
#' fpca.preprocess converts discrete point observations of a functional time series into functional time series evaluated on a dense
#' workinggrid. It also computes its eigendecomposition and empirical Karhunen-Loeve expansion. The data is preprocessed using
#' either natural splines, smoothing splines, or the PACE-method.
#'
#' @param data A multivariate time series of class 'ts' or 'data.frame'.
#' @param workinggrid An optional list, sequence or numeric vector which contains a predetermined workinggrid for the analysis.
#' If NULL, the sequence with gridsize 1 is used, ranging from the smallest to the highest observed domain point. NULL is default.
#' @param observationgrid An optional list, sequence or numeric vector which contains a predetermined observationgrid for the
#' observations. If NULL, the column names will be used. NULL is default.
#' @param method Preprocessing method to be used. If 'smoothsplines', smoothing splines from the package 'fda' are used (see Kokoszka
#' and Reimherr, 2017). If 'pace', the PACE method from the package 'fdapace' is applied (see Yao et al., 2005).
#' If 'natralsplines', natural spline interpolation is used.
#' In case of irregularly observed functional data please set to 'pace' or 'naturalsplines'.
#' 'smoothsplines' is default.
#'
#' @return
#' fpca.preprocess returns an object of class 'FPCAobj' with:
#' \item{eigenvalues}{A vector containing the eigenvalues of the sample covariance operator.}
#' \item{fpcascores}{A matrix/ts-matrix of the computed FPCA-scores/projection coefficients.}
#' \item{eigenfunctions.obsgrid}{A matrix containing the orthonormal eigenfunctions of the sample covariance operator, which are evaluated on the workinggrid.}
#' \item{meanfunction.obsgrid}{A vector containing the sample mean function, which is evalauted on the observationgrid.}
#' \item{eigenfunctions.workgrid}{A matrix containing the orthonormal eigenfunctions of the sample covariance operator, which are evaluated on the observationgrid.}
#' \item{meanfunction.workgrid}{A vector containing the sample mean function, which is evalauted on the workinggrid.}
#' \item{raw.data}{The raw data used for the analysis.}
#' \item{observationgrid}{Observationgrid of the raw data used for the analysis.}
#' \item{workinggrid}{The dense workinggrid of the analysis.}
#' \item{dense.fts}{The dense representation of the functional time series using the empirical Karhunen-Loeve expansion on the workinggrid.}
#' \item{method}{Preprocessing method to obtain the dense functional time series.}
#'
#' @export
#' @references
#' * Kokoszka P. and Reimherr M. (2017), "Introduction to Functional Data Analysis". CRC Press.
#' * Otto S. and Salish N. (2022). "Dynamic Factor Model for Functional Time Series: Identification, Estimation, and Prediction". ...
#' * Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). Functional data analysis for sparse longitudinal
#' data. Journal of the American Statistical Association, 100:577-590.
#' @examples
#' # standard workinggrid
#' fed = load.fed()
#' fpca.preprocess(data = fed)
#' # A tighter workinggrid
#' wg = (2:720)/2
#' fpca.preprocess(fed, workinggrid = wg)
fpca.preprocess = function(data, workinggrid = NULL, observationgrid = NULL, method = c("smoothsplines", "pace", "naturalsplines")){
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
  if(!all(diff(observationgrid)>0))stop("observationgrid must hyve strictly increasing entries")
  observationgrid.unit = (observationgrid-head(observationgrid,1))/(tail(observationgrid,1) - head(observationgrid,1))
  # define equidistant workinggrid based on minimum difference in the observationgrid
  if(is.null(workinggrid)) workinggrid = seq(observationgrid[1], observationgrid[length(observationgrid)], min(diff(observationgrid)))
  workinggrid.unit = (workinggrid-head(workinggrid,1))/(tail(workinggrid,1) - head(workinggrid,1))
  ## analysis
  # naturalsplines
  if(method == "naturalsplines"){
    nsfit = ts(matrix(nrow = dim(data)[1], ncol = length(workinggrid)), start=time(data)[1], frequency = frequency(data))
    colnames(nsfit)=workinggrid
    for(i in 1:dim(data)[1]){
      thisobsgrid = observationgrid[!is.na(data[i,])]
      nsbasis = splines::ns(thisobsgrid, knots=thisobsgrid[-length(thisobsgrid)], Boundary.knots=c(thisobsgrid[1],thisobsgrid[length(thisobsgrid)]))
      coef = lm(na.omit(as.matrix(data)[i,]) ~ nsbasis - 1)$coefficients
      densebasis = splines::ns(workinggrid, knots=thisobsgrid[-length(thisobsgrid)], Boundary.knots=c(thisobsgrid[1],thisobsgrid[length(thisobsgrid)]))
      nsfit[i,] = densebasis %*% coef
    }
    ## FPCA
    pca = prcomp(nsfit)
    ## Norms of eigenfunctions are approximated on the workinggrid using the trapezoidal rule for numerical intergation
    ## eigenfunctions are normalized and corresponding scores defined accordingly
    L2norm = function(z) pracma::trapz(workinggrid, z^2)
    norm.eigenf = apply(pca$rotation, 2, L2norm)
    sig.eigenf=sign(pca$rotation[1,])
    sig.eigenf[which(sig.eigenf==0)]=1
    eigenfunctions.workgrid = t(c(sig.eigenf/sqrt(norm.eigenf)) * t(pca$rotation))
    fpcascores = t(c(sig.eigenf*sqrt(norm.eigenf)) * t(pca$x))
    if(is.ts(data)) fpcascores = ts(fpcascores, start = head(time(data),1), frequency = frequency(data))
    eigenvalues = norm.eigenf*pca$sdev^2
    meanfunction.workgrid = pca$center
    ##
    eigenfunctions.obsgrid = eigenfunctions.workgrid[match(observationgrid, workinggrid),]
    meanfunction.obsgrid = meanfunction.workgrid[match(observationgrid, workinggrid)]
  }
  # smoothingsplines
  if(method == "smoothsplines"){
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

  dense.fts = fpcascores %*% t(eigenfunctions.workgrid) + matrix(rep(meanfunction.workgrid, dim(fpcascores)[1]), nrow = dim(fpcascores)[1], byrow=TRUE)
  if(is.ts(data)) dense.fts = ts(dense.fts, start = head(time(data),1), frequency = frequency(data))

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
    "dense.fts" = dense.fts,
    "method" = method
  )
  class(out) = "FPCAobj"
  return(out)
}
