#' load Liu-Wu Yield Data
#'
#' Function which provides data from Liu-Wu Yield Data via the following website:
#' https://sites.google.com/view/jingcynthiawu/yield-data
#' 
#' Yan Liu and Jing Cynthia Wu "Reconstructing the Yield Curve", Journal of Financial Economics, 2021, 142 (3), 1395-1425. 
#'
#' @return
#' Annualized continuously-compounded zero-coupon yields in percentage points. Each column corresponds to a maturity between 1 to 360 months.
#' @export
#' @examples
#' data = load.LW()
load.LW <- function(){
  data = gsheet::gsheet2tbl('https://docs.google.com/spreadsheets/d/15uGZet8bS5rEHy9nbGmgCHXyMGUSoHx3/edit?usp=drive_link&ouid=117915996921355706819&rtpof=true&sd=true')
  df = data[-(1:8), -1]
  for(i in seq_along(df)) {
    if(is.character(df[[i]])) {
      df[[i]] <- as.numeric(as.character(df[[i]]))
    }
  }
  LW.data = ts(df, start = c(1961,6), frequency=12)
  colnames(LW.data) = 1:360
  LW.data = window(LW.data, start = 1986, end = c(2023,12))
  return(LW.data)
}