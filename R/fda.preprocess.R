#' @import stats
#' @import utils
NULL

#' Natural Spline Preprocessing of Funtional Data
#'
#' The `fda.preprocess` function transforms discrete point observations of functional data (functional time series) into functional data evaluated on a dense equidistant grid. This function uses a direct interpolation approach, using natural splines to interpolate missing values. In addition, it performs the eigendecomposition of the sample covariance operator and computes the empirical Karhunen-Loève expansion (empirical functional principal components). Missing values at the beginning or end are reconstructed using the Kneip and Liebl (2020) optimal reconstruction operator.
#'
#' @param data A multivariate time series or data in 'ts' or 'data.frame' format. The function expects data to represent discrete observations of functional data (functional time series) that may contain missing values.
#' @param observationgrid An optional numeric vector specifying the observation grid of the input data. If NULL (default), the function attempts to infer the grid from column names of 'data'.
#' @param workinggrid An optional numeric vector defining an equidistant grid for the analysis. If NULL (default), the function generates an equidistant grid based on the minimum difference in 'observationgrid'
#'
#' @return
#' Returns an object of class 'fdaobj', which includes the following components:
#' \item{densedata}{The data interpolated onto the dense working grid.}
#' \item{workinggrid}{The generated or specified dense equidistant working grid.}
#' \item{operator}{The operator for which the eigenelements are computed. Here: "sample_covariance(FPC)".}
#' \item{scores}{The coefficients representing the projections of 'densedata' onto each eigenfunction.}
#' \item{eigenfunctions}{A matrix of orthonormal eigenfunctions derived from the sample covariance operator, representing the functional principal components.}
#' \item{eigenvalues}{The eigenvalues associated with each eigenfunction, indicating the variance explained by each principal component.}
#' \item{scores.centered}{The projection coefficients after demeaning 'densedata'.}
#' \item{meanfunction}{The sample mean function calculated across the dense data.}
#' \item{raw.data}{The original input data.}
#' \item{observationgrid}{The observation grid used or inferred from the input data.}
#'
#' @export
#' @references
#' * Hsing, T., & Eubank, R. (2015). Theoretical foundations of functional data analysis, with an introduction to linear operators. John Wiley & Sons.
#' * Kneip, A., & Liebl, D. (2020). On the optimal reconstruction of partially observed functional data. The Annals of Statistics, 48, 1692-1717.
#' * Otto, S., & Salish, N. (2024). Approximate Factor Models For Functional Time Series. arXiv:2201.02532.
#' @examples
#' # Example with standard working grid
#' fed = load.fed()
#' fda.preprocess(data = fed)
#'
#' # Example with a customized tighter working grid
#' wg = seq(1,360, by=0.5)
#' fda.preprocess(fed, workinggrid = wg)
fda.preprocess = function(data, observationgrid = NULL, workinggrid = NULL){
  ## data: A multivariate time series of class 'ts' or 'data.frame'.
  ## workinggrid: An optional vector containing an equidistant workinggrid for the analysis. If NULL, an equidistant workinggrid is used with gridsize equal to the minimum difference in the observationgrid. NULL is default.
  ## observationgrid: An optional vector which contains a predetermined observationgrid for the observations. If NULL, the column names will be used. NULL is default.
  ## #############################################################################
  ## check if input is valid and define observationgrid and workinggrid if needed
  ## #############################################################################
  ## check data
  if(missing(data) || !is.numeric(data) || !is.matrix(data)){
    stop("Please provide valid data. See 'data' in the help file for more details.")
  }
  ## check observationgrid
  if(is.null(observationgrid)) observationgrid = as.numeric(colnames(data))
  if(!is.vector(observationgrid)) stop('observationgrid must be a vector. Please specify a valid observation grid')
  if(is.vector(observationgrid)){
    if(length(observationgrid) != ncol(data)) stop('Length of observationgrid does not fit. Please specify a valid observation grid.')
    if(!is.numeric(observationgrid)) stop("observationgrid must be numeric. Please specify a valid observation grid.")
    if(any(is.na(observationgrid))) stop('observationgrid contains NA. Please specify a valid observation grid.')
    if(!all(diff(observationgrid)>0)) stop("observationgrid must be strictly increasing. Please specify a valid observation grid.")
  }
  ## check workinggrid
  if(!is.null(workinggrid)){
    work.diff <- diff(workinggrid)
    if(all(abs(work.diff - work.diff[1]) > 1e-8)){
      stop("workinggrid must be equidistant. Please specify a valid working grid.")
    }
    if(!all(work.diff>0)){
      stop("workinggrid must be strictly increasing. Please specify a valid working grid.")
    }
  } else {
    ## define equidistant workinggrid based on minimum difference in the observationgrid
    workinggrid = seq(observationgrid[1], observationgrid[length(observationgrid)], min(diff(observationgrid)))
  }

  ## #############################################################################
  ## preprocess using natural splines and optimal reconstruction
  ## #############################################################################
  ## Fit natural splines to inperpolate missing data
  splinefit = matrix(nrow = dim(data)[1], ncol = length(workinggrid), dimnames = list(NULL,workinggrid))
  for(i in 1:dim(data)[1]){
    ## define a natural spline basis with knots at the observed grid of curve i
    thisobsgrid = observationgrid[!is.na(data[i,])]
    splinebasis = splines::ns(thisobsgrid, knots=thisobsgrid[-length(thisobsgrid)], Boundary.knots=c(thisobsgrid[1],thisobsgrid[length(thisobsgrid)]))
    ## fit natural splines
    coef = lm(na.omit(as.matrix(data)[i,]) ~ splinebasis - 1)$coefficients
    ## compute spline fit on dense workinggrid restricted to the observed range of curve i
    thisrange = workinggrid >= thisobsgrid[1] & workinggrid <= rev(thisobsgrid)[1]
    densebasis = splines::ns(workinggrid[thisrange], knots=thisobsgrid[-length(thisobsgrid)], Boundary.knots=c(thisobsgrid[1],rev(thisobsgrid)[1]))
    splinefit[i,thisrange] = densebasis %*% coef
  }
  ## Reconstruct missing parts at the beginning/end of domain
  if(any(is.na(splinefit))){
    ffm.reconstruct = function(dat){
      Ly = lapply(seq_len(nrow(dat)), function(i) dat[i,])
      Lu = lapply(seq_len(nrow(dat)), function(i) as.numeric(colnames(dat)))
      reconst_result = ReconstPoFD::reconstructKneipLiebl(Ly = Ly, Lu = Lu,
                                                          method = 'Error=0_AlignYES_CommonGrid', reconst_fcts = NULL)
      dat.reconst = matrix(unlist(reconst_result[['Y_reconst_list']]), ncol=dim(dat)[2], byrow=TRUE)
    }
    message("Your curve data contains missing parts at the beginning or end of the domain that need to be reconstructed. The 'ReconstPoFD' method from Kneip and Liebl (2020) is applied. This might take a while.")
    if(!requireNamespace("ReconstPoFD")){
      stop("Please install 'ReconstPoFD' from https://github.com/lidom/ReconstPoFD and run your code again..")
    }
    splinefit = ffm.reconstruct(splinefit)
  }
  if(is.ts(data)){
    splinefit = ts(splinefit, start = time(data)[1], frequency=frequency(data))
  }
  dimnames(splinefit) = list(1:dim(splinefit)[1], round(workinggrid,2))

  ## #############################################################################
  ## compute empirical Karhunen-Loeve decomposition
  ## #############################################################################

  FPCA = fda.FPCA(splinefit, workinggrid)

  ## output
  out=list(
    "densedata" = splinefit,
    "workinggrid" = workinggrid,
    "operator" = "sample_covariance(FPC)",
    "scores" = FPCA$scores,
    "eigenfunctions" = FPCA$eigenfunctions,
    "eigenvalues" = FPCA$eigenvalues,
    "scores.centered" = FPCA$scores.centered,
    "meanfunction" = FPCA$meanfunction,
    "raw.data" = data,
    "observationgrid" = observationgrid
  )
  class(out) = "fdaobj"
  return(out)
}




fda.FPCA = function(densedata, workinggrid, start = NULL, end = NULL){
  ## ####################################################################
  ## Compute KL decomposition (eigenelements of the covariance operator)
  ## ####################################################################
  ## workinggrid must be equidistant;
  ## densedata must be fully observed on workinggrid;
  ## start and end are optional start and end indices and define the
  ## restricted period for which the eigendecomposition is computed.
  ## ####################################################################
  if(dim(densedata)[2] != length(workinggrid)){
    stop("Please specify a valid workinggrid.")
  }
  ## check start and end, restrict factors to specified period
  if(is.null(start)){
    start = 1
  } else {
    if((start < 1) || (start > dim(densedata)[1])){
      stop("Please specify valid start.")
    }
  }
  if(is.null(end)){
    end = dim(densedata)[1]
  } else {
    if((end < 1) || (end > dim(densedata)[1])){
      stop("Please specify valid end")
    }
  }
  ## restrict data to specified period
  if(is.ts(densedata)){
    redata = window(densedata, start = time(densedata)[start], end = time(densedata)[end])
  } else {
    redata = window(densedata, start = start, end = end)
  }
  ## sample covariance matrix for redata
  n=dim(redata)[1]
  C = cov(redata)*(n-1)/n
  eig = eigen(C)
  ## signs of first entry of eigenvectors
  signs = ifelse((sign(eig$vectors[1,])==-1), -1, 1)
  ## binwidth for rectangular numerical integration based on workinggrid
  binwidth = (rev(workinggrid)[1]-(workinggrid[1]))/length(workinggrid)
  ## flip signs and normalize so that L2norm is 1
  eigenfunctions = t(signs*t(eig$vectors))/sqrt(binwidth)
  ## eig$values corresponds to binwidth=1 (matrixproduct)
  ## multiply with given binwidth to get eigenvalues in L2-sense
  eigenvalues = eig$values*binwidth
  ## compute projection coefficients (scores) for full densedata
  scores = (densedata %*% eigenfunctions)*binwidth
  ## centered scores and meanfunction (KL)
  scores.centered = scores - matrix(rep(colMeans(scores), dim(densedata)[1]), nrow = dim(densedata)[1], byrow=TRUE)
  meanfunction = colMeans(densedata)

  ## #############################################################################
  ## define column names and transform back to ts if input was ts
  ## #############################################################################
  dimnames(eigenfunctions) = list(workinggrid, paste0("FPC", 1:length(workinggrid)))
  dimnames(scores) = list(1:dim(densedata)[1], paste0("FPC", 1:length(workinggrid)))
  dimnames(scores.centered) = list(1:dim(densedata)[1], paste0("FPC", 1:length(workinggrid)))
  if(is.ts(densedata)){
    scores = ts(scores, start = time(densedata)[1], frequency=frequency(densedata))
    scores.centered = ts(scores.centered, start = time(densedata)[1], frequency=frequency(densedata))
  }
  out=list(
    "densedata" = densedata,
    "scores" = scores,
    "eigenfunctions" = eigenfunctions,
    "eigenvalues" = eigenvalues,
    "scores.centered" = scores.centered,
    "meanfunction" = meanfunction,
    "workinggrid" = workinggrid
  )
  return(out)
}

