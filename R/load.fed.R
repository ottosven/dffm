#' load fed
#'
#' Function which provides data from the FED via the following website:
#' https://www.federalreserve.gov/datadownload/Output.aspx?rel=H15&series=bf17364827e38702b42a58cf8eaa3f78&lastobs=&from=&to=&filetype=csv&label=include&layout=seriescolumn&type=package
#'
#' @return
#' Time series of yield curve data containing 11 maturities ranging from 1 to 360 month, which are
#' observed from 2002 to 2022.
#' @export
#' @examples
#' data = load.fed()
load.fed <- function(){
  temp <- tempfile()
  download.file("https://www.federalreserve.gov/datadownload/Output.aspx?rel=H15&series=bf17364827e38702b42a58cf8eaa3f78&lastobs=&from=&to=&filetype=csv&label=include&layout=seriescolumn&type=package",temp)
  data = read.csv(temp, skip = 5)
  unlink(temp)
  dates = data[,1]
  options(warn = -1)
  data[,2] = as.numeric(data[,2])
  data[,3] = as.numeric(data[,3])
  data[,4] = as.numeric(data[,4])
  data[,5] = as.numeric(data[,5])
  data[,6] = as.numeric(data[,6])
  data[,7] = as.numeric(data[,7])
  data[,8] = as.numeric(data[,8])
  data[,9] = as.numeric(data[,9])
  data[,10] = as.numeric(data[,10])
  data[,11] = as.numeric(data[,11])
  data[,12] = as.numeric(data[,12])
  options(warn = 0)
  colnames(data) = c("dates",1,3,6,12,24,36,60,84,120,240,360)
  completeVec = complete.cases(data)
  data = data[completeVec,]
  data = ts(data[,-1], frequency = 12, start = c(2001,7), end=c(2021,12))
  return(suppressWarnings(data))
}
