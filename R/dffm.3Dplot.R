#' 3D plot of Functional Time Series
#'
#' A function that produces a 3D-plot of the data used in the functions 'fda.preprocess'.
#'
#' @param fpcaobj Must be an object of class 'fdaobj'.
#' @param domainlab ylab name of the plot.
#' @param outputlab zlab name of the plot.
#' @param rotate Parameter to change the horizontal rotation of the 3D plot. If NULL rotate will be set to 35.
#' NULL is default.
#' @param cex.main Parameter to change the font size of the headline. If NULL cex.main will be set to 1.6.
#' NULL is default.
#' @param ... Option to change parameters in surf3D().
#'
#' @return
#' 3D plot of data used in the functions 'fda.preprocess'.
#' @export
#'
#' @examples
#' JKV = load.LW()
#' d = fda.preprocess(JKV)
#' d$raw.data
#' dffm.3Dplot(d)
dffm.3Dplot = function(fpcaobj, domainlab=NULL, outputlab=NULL, rotate = NULL, cex.main = NULL, ...){
  if(is.null(rotate)) rotate = 35
  if(is.null(cex.main)) cex.main = 1.6
  observationgrid = fpcaobj$observationgrid
  data = fpcaobj$raw.data
  Tdim = dim(data)[1]
  time = replicate(length(observationgrid), c(time(data)))
  plot3D::surf3D(x = time,
                 y = t(replicate(Tdim, observationgrid)),
                 z = data,
                 bty="b2",
                 theta=rotate,
                 phi=12,
                 colvar = data,
                 col = NULL,
                 shade = 0.05,
                 colkey = FALSE,
                 border = "black",
                 zlim = range(data) + c(-0.1,0.1),
                 xlab = 'Time',
                 ylab = domainlab,
                 zlab = outputlab,
                 ticktype = "detailed",
                 cex.main=cex.main,
                 expand = 0.6,
                 ...
  )
}
