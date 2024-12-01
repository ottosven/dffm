#' load jungbacker
#'
#' A function which takes yield curve data provided by Jungbacker et.al (2014) from the
#' following website: http://qed.econ.queensu.ca/jae/2014-v29.1/jungbacker-koopman-van_der_wel/
#'
#' @return
#' Time series of yield curve data containing maturities ranging from 3 to 120 month, which had
#' been observed from 1987 to 2008.
#' @references
#' * Jungbacker, B., Koopman, S. J., and Van der Wel, M. (2014). Smooth dynamic factor analysis with application to the US term structure of
#' interest rates. Journal of Applied Econometrics, 29:65-90.
#' @export
#' @examples
#' data = load.JKV()
load.JKV = function(){
  temp = tempfile()
  download.file("http://qed.econ.queensu.ca/jae/2014-v29.1/jungbacker-koopman-van_der_wel/jkv-data.zip",temp)
  data = read.table(unz(temp, "UnsmFB_70-09.txt"))
  unlink(temp)
  jungbacker.data = ts(data, start = c(1970, 1), frequency = 12)
  colnames(jungbacker.data) = c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
  jungbacker.data =  window(jungbacker.data, start=c(1987,1), end=c(2007,12))
  return(jungbacker.data)
}
