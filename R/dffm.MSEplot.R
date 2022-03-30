#' dffm.MSEplot
#'
#' A function which plots mean squared errors from the dffm.criterion in 3D. For more information about the
#' mse errors, see the dffm.criterion function.
#'
#' @param criterion Has to be an object of class 'dffm', 'FPCAobj' or results of the dffm.criterion function.
#' @param kmax A parameter specifying the factors used in the 3D plot.
#' @param pmax A parameter specifying the lags used in the 3D plot.
#' @param ... Option to change parameters in surf3D().
#'
#' @return
#' A 3D plot of the mse errors from the dffm.criterion function.
#' @export
#'
#' @examples
#' # with dffm.criterion
#' JKV = load.JKV()
#' fpca = fpca.preprocess(data = JKV, method = "splines")
#' dffm.MSEplot(dffm.criterion(fpca, p.max = 8), kmax = 15, pmax = 8)
#' # with FPCAobj
#' fpca = fpca.preprocess(data = JKV, method = "splines")
#' dffm.MSEplot(fpca, kmax = 12, pmax = 10)
dffm.MSEplot = function(criterion, kmax, pmax, ...){
  if(!"MSE.matrix" %in% names(criterion) & class(criterion) == "dffm" | !"MSE.matrix" %in% names(criterion) & class(criterion) == "FPCAobj"){
    criterion = dffm.criterion(criterion)
  }
  if(!"MSE.matrix" %in% names(criterion) & class(criterion) == "dffm" & class(criterion) == "FPCAobj"){
    stop("criterion has to be of class 'dffm' or 'FPCAobj' or an output of the function dffm.criterion")
  }
  mse = criterion$MSE.matrix[1:kmax, 1:pmax]
  plot3D::surf3D(y = replicate(pmax, 1:kmax),
                 x = t(replicate(kmax, 1:pmax)),
                 ytick<-seq(1, 10, by=5),
                 z = mse,
                 bty="b2",
                 theta=120,
                 phi=12,
                 colvar = mse,
                 col = plot3D::jet2.col(n=10^3),
                 shade = 0.05,
                 colkey = FALSE,
                 border = "black",
                 xlab = 'Number of lags',
                 ylab ='Number of factors',
                 zlab ='MSE',
                 ticktype = "detailed",
                 cex.main=1.2,
                 cex.lab=1,
                 cex.axis=0.8,
                 expand = 0.6,
                 ...
  )
}
