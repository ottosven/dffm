#' Plot of the MSE of different combinations of factors and lags
#'
#' A function that produces a 3D-plot of the mean squared errors from the dffm.criterion function. For more information about the
#' mse, see the fts.criterion function.
#'
#' @param criterion Must be an output of the fts.criterion function.
#' @param kmax A parameter specifying the maximum number of factors used in the 3D plot.
#' @param pmax A parameter specifying the maximum number lags used in the 3D plot.
#' @param rotate Parameter to change the horizontal rotation of the 3D plot. If NULL rotate will be set to 120.
#' NULL is default.
#' @param cex.main Parameter to change the font size of the headline. If NULL cex.main will be set to 1.2.
#' NULL is default.
#' @param ... Option to change parameters in surf3D().
#'
#' @return
#' A 3D plot of the mse errors from the dffm.criterion function.
#' @export
#'
#' @examples
#' # with dffm.criterion
#' JKV = load.JKV()
#' JKV.fda = fda.preprocess(data = JKV)
#' dffm.MSEplot(fts.criterion(JKV.fda, p.max = 8), kmax = 15, pmax = 8)
dffm.MSEplot = function(criterion, kmax, pmax, rotate = NULL, cex.main = NULL, ...){
  if(sum(which(criterion$MSE.matrix[1:kmax, 1:pmax] == "Inf", arr.ind = TRUE, useNames = TRUE)) > 0){
    stop("criterion is not allowed to have 'Inf' values inside of the 'MSE.matrix'")
  }
  if(is.null(rotate)) rotate = 120
  if(is.null(cex.main)) cex.main = 1.2
  mse = criterion$MSE.matrix[1:kmax, 1:pmax]
  plot3D::surf3D(y = replicate(pmax, 1:kmax),
                 x = t(replicate(kmax, 1:pmax)),
                 ytick<-seq(1, 10, by=5),
                 z = mse,
                 bty="b2",
                 theta=rotate,
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
                 cex.main=cex.main,
                 expand = 0.6,
                 ...
  )
}
