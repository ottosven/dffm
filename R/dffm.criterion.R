#' dffm.criterion
#'
#' Estimates the number of factors K and lags p for the functional factor model using the
#' information criteria proposed in Otto and Salish (2021).
#'
#' @param fpcaobj An object of class 'FPCAobj' or 'dffm'.
#' @param K.max The maximum number of factors for the criterion. If NULL K.max will be highest possible number from the given 'FPCAobj'.
#' NULL is default.
#' @param p.max The maximum number of lags for the criterion. If NULL p.max will be based on Schwert's rule of thumb. NULL is default.
#'
#' @return
#' A list containing the following:
#' \item{IC.min}{Estimated number of factors K and lags p, based on the BIC, HQC and fFPE criterion.}
#' \item{MSE.matrix}{A matrix containing all computed MSE combinations.}
#'
#' @export
#' @references
#' * Otto S. and Salish N. (2022). "Dynamic Factor Model for Functional Time Series: Identification, Estimation, and Prediction". ...
#'
#' @examples
#' # using object of class 'FPCAobj'
#' fed = load.fed()
#' fpca = fpca.preprocess(data = fed)
#' dffm.criterion(fpca)
#' # using object of class 'dffm'
#' d = dffm(fed)
#' dffm.criterion(d)
dffm.criterion = function(fpcaobj, K.max = NULL, p.max = NULL){
  # checking data
  if(class(fpcaobj)[1]!="FPCAobj" & class(fpcaobj)[1] == "dffm"){
    fpcaobj = fpcaobj$FPCA
  }
  if(class(fpcaobj)[1] != "FPCAobj" & class(fpcaobj)[1] != "dffm")stop("'fpcaobj' has to be an object of class 'FPCAobj' or 'dffm'")
  N = dim(fpcaobj$fpcascores)[1]
  if(is.null(K.max)) (K.max = dim(fpcaobj$fpcascores)[2])
  if(is.null(p.max)) (p.max = max(floor(8*(N/100)^(1/4)),1)) #Schwert rule of thumb
  if(K.max > length(fpcaobj$eigenvalues)) stop(paste("K.max is not allowed to be greater then", length(fpcaobj$eigenvalues)))
  ## Evaluate IC
  eval.criterion = function(J,m){
    if(J*m >= (N-J*m)) (return(rep(Inf,4)))
    factors = fpcaobj$fpcascores[,1:J, drop=FALSE]
    colnames(factors) = paste("F.",1:J, sep = "")
    VARfit = lm(embed(factors, m+1)[,1:J] ~ -1 + embed(factors, m+1)[,-(1:J)])
    if(VARfit$rank != J*m) (return(rep(Inf,4)))
    tr.sigmahat = sum((residuals(VARfit))^2)/(N-(J*m))
    tot.prederror = tr.sigmahat + sum(fpcaobj$eigenvalues[-c(1:J)])
    CR = c(
      tot.prederror,
      log(tot.prederror) + J*m*log(N)/N,
      log(tot.prederror) + 2*J*m*log(log(N))/N,
      (N+m*J)/N * tr.sigmahat + sum(fpcaobj$eigenvalues[-c(1:J)])
    )
    return(CR)
  }
  IC=list()
  IC.min=matrix(nrow = 3, ncol=2, dimnames = list(c("BIC", "HQC", "fFPE"), c("K.opt", "p.opt")))
  for(j in 1:4) IC[[j]] = matrix(nrow = K.max, ncol = p.max, dimnames = list(paste("K.", 1:K.max, sep = ""), paste("p.", 1:p.max, sep = "")))
  for (i in 1:K.max){
    aux = sapply(1:p.max, eval.criterion, J = i)
    for(j in 1:4) IC[[j]][i,] = aux[j,]
  }
  for(j in 1:3) IC.min[j,] = which(IC[[j+1]] == min(na.omit(IC[[j+1]])), arr.ind = TRUE)[1,]
  MSE.matrix = IC[[1]]
  return(list("IC.min" = IC.min,"MSE.matrix" = MSE.matrix))
}
