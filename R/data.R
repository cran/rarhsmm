#' Pseudo activity count (per minute) data for cats
#'
#' A dataset containing the daily closing price of 50 NYSE stocks
#' from 2014-01-10 to 2017-01-01
#'
#' @format A data frame with 756 rows and 50 variables:
## \describe{
##   \item{id}{cat ID: 1,2,3}
##   \item{hour}{hour of the day: 1,2,...,24}
##   \item{minute}{minute of the hour: 1,2,...,60}
##   \item{night}{night time indicator}
##   \item{activity}{activity count data}
## }
"finance"

#library(tseries)
#?get.hist.quote
#stocks <- c("AMD","BAC","AAPL","FTR","BBRY","GE","RAD","MU","CHK","F",
#            "VALE","FCX","PBR","XOM","INTC","MSFT","WWAV","QQQ","ABEV",
#            "VZ","HPQ","KKFC","PFE","SWN","T","AKS","FCAU","JPM","SFUN",
#            "LULU","MT","WFT","CLF","SNA","S","C","HCP","DRYS","FMC",
#            "CSCO","KMI","AES","X","SIRI","WLL","COP","BSX","JCP","RF","EDUC")

#list1 <- vector(mode = "list", length = 50)

#for(i in 1:50){
#  list1[[i]] <- get.hist.quote(instrument=stocks[i], start="2014-01-01",
#                               end="2017-01-01",quote="Close")
#}

#sapply(list1,length)

#finance <- matrix(0,length(list1[[1]]),50)
#for(i in 1:50) finance[,i] <- as.vector(list1[[i]])
#finance <- data.frame(finance)
#names(finance) <- stocks

#take log returns!!!!!

#save(finance,file="stocks_150101_170101.RData")


