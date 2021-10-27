#' @import stats
#' @import fda
#' @import xts
#' @import tsbox
#' @import fdapace
#' @import vars
#' @import zoo
NULL

#' predict dffm
#'
#' The generic S3 method predict of an object of class 'dffm' will make predictions based on vector autoregressive
#' or autoregressive models.
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
#' predict(d)
predict.dffm = function(x, AR = FALSE, criterion = FALSE, h = NULL, p = NULL, K = NULL){
  if(class(x) != "dffm") stop('x must be class "dffm"')
  if(criterion == TRUE & !is.null(p) | criterion == TRUE & !is.null(K)| criterion == TRUE & !is.null(K) & !is.null(p)){
    stop("please use either criterion or predetermined K and/or p")}
  if(criterion == TRUE){
    IC = dffm.criterion(x)
    p = as.numeric(IC$IC.min[1])
    K = as.numeric(IC$IC.min[2])
  }
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
  colnames(fitt) = x$basis$workgrid
  fitted.values = fitt
  output = list("predicted.values" = fitted.values, "predicted.factors" = f.cast)
  return(output)
}
