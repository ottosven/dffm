library(dffm)
library(demography)
library(gsheet)


# Create data directory if it doesn't exist
dir.create("data", showWarnings = FALSE)

##  ###########################################################################
##  Mortality data
##  Registration on the website https://www.mortality.org is required.
username = ""
password = ""
##  ###########################################################################


USA.raw = hmd.mx("USA", username, password, "USA")
USA = ts(log(t(USA.raw$rate$male)), start = 1933, frequency=1)
saveRDS(USA, file = "./data/mort.rds")
saveRDS(fda.preprocess(USA), file = "./data/mort-fda.rds")


##  ###########################################################################
##  Yield curve data
##  Available on Jing Cynthia Wu's website: https://sites.google.com/view/jingcynthiawu/yield-data
##  ###########################################################################

data = gsheet::gsheet2tbl('https://docs.google.com/spreadsheets/d/1-wmStGZHLx55dSYi3gQK2vb3F8dMw_Nb/edit?gid=1606479329#gid=1606479329')
df = data[-(1:8), -1]
for(i in seq_along(df)) {
  if(is.character(df[[i]])) {
    df[[i]] = as.numeric(as.character(df[[i]]))
  }
}
LW.data = ts(df, start = c(1961,6), frequency=12)
colnames(LW.data) = 1:360
LW360 = window(LW360, start = c(1985,11), end = c(2023,12))
saveRDS(LW360, "./data/LW360.rds")
saveRDS(fda.preprocess(LW360), "./data/LW360-fda.rds")
