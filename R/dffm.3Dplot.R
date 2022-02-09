#' dffm.3Dplot
#'
#' A function which 3D-plots the data used in the functions 'fpca.preprocess' or 'dffm'. For more details
#' look up fpca.preprocess() or dffm().
#'
#' @param fpcaobj Has to be an object of class 'FPCAobj' or 'dffm'.
#' @param domainlab ...
#' @param outputlab ...
#' @param ... Option to change parameters in surf3D().
#'
#' @return
#' 3D plot of data used in the functions 'fpca.preprocess' or 'dffm'.
#' @export
#'
#' @examples
#' JKV = load.JKV()
#' d = dffm(JKV)
#' d$FPCA$raw.data
#' dffm.3Dplot(d)
dffm.3Dplot = function(fpcaobj, domainlab=NULL, outputlab=NULL, ...){
  if(class(fpcaobj) == "dffm") fpcaobj = fpcaobj$FPCA
  if(class(fpcaobj) != "dffm" & class(fpcaobj) != "FPCAobj")stop("fpcaobj has to be an object of class 'FPCAobj' or 'dffm'")
  observationgrid = fpcaobj$observationgrid
  data = fpcaobj$raw.data
  Tdim <- dim(data)[1]
  time <- replicate(length(observationgrid), c(time(data)))
  plot3D::surf3D(x = time,
                 y = t(replicate(Tdim, observationgrid)),
                 z = data,
                 bty="b2",
                 theta=35,
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
                 cex.main=1.6,
                 cex.lab=1,
                 cex.axis=1,
                 expand = 0.6,
                 ...
  )
}
