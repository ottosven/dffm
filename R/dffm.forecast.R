#' dffm.forecast
#'
#' A function to predict h month ahead curves from a 'dffm' object.
#'
#' @param dffmobj An object of class 'dffm'.
#' @param h A parameter determining the forecast horizon. If h is NULL, h will be 1. NULL is default.
#'
#' @return
#' A list containing the following:
#' \item{predfactors}{A vector containing the predicted factors for further computations.}
#' \item{predcurves.obsgrid}{A matrix or ts-matrix containing the predicted curves evaluated on the observationgrid.}
#' \item{predcurves.workgrid}{A matrix or ts-matrix containing the predicted curves evaluated on the workinggrid.}
#' @export
#'
#' @examples
#' fed = load.fed()
#' d = dffm(fed, K = 3, p = 2, AR = FALSE)
#' dffm.forecast(d, h = 2)
dffm.forecast = function(dffmobj, h=1){
  if(class(dffmobj) != "dffm")stop("'dffmobj' has to be an object of class 'dffm'")
  allfactors = dffmobj$factors
  for(i in 1:h) allfactors = rbind(allfactors, tail(embed(allfactors,dffmobj$p),1) %*% dffmobj$VARcoefficients)
  if(is.ts(dffmobj$factors)) allfactors = ts(allfactors, start = start(dffmobj$factors), frequency = frequency(dffmobj$factors))
  predfactors = window(allfactors,start = tsp(dffmobj$factors)[2]+1/(frequency(dffmobj$factors)))
  predcurves.obsgrid = predfactors %*%  t(dffmobj$loadingfunctions.obsgrid) + matrix(rep(t(dffmobj$meanfunction.obsgrid),h), nrow = h, byrow=TRUE)
  predcurves.workgrid = predfactors %*%  t(dffmobj$loadingfunctions.workgrid) + matrix(rep(t(dffmobj$meanfunction.workgrid),h), nrow = h, byrow=TRUE)
  if(is.ts(dffmobj$factors)){
    predcurves.obsgrid = ts(predcurves.obsgrid, start=start(predfactors), frequency =  frequency(predfactors))
    predcurves.workgrid = ts(predcurves.workgrid, start=start(predfactors), frequency =  frequency(predfactors))
  }
  out=list(
    "predfactors" = predfactors,
    "predcurves.obsgrid" = predcurves.obsgrid,
    "predcurves.workgrid" = predcurves.workgrid
  )
  return(out)
}
