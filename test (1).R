#test data

library(quantmod) # use for gathering and charting economic data
library(lubridate) # date functions

par(mfrow = c(2,2)) # four plots on one window/page

# Economic Data from Federal Reserve Bank of St. Louis (FRED system)
# National Civilian Unemployment Rate (monthly, percentage)
getSymbols("UNRATENSA", src="FRED", return.class = "xts")
ER <- 100 - UNRATENSA # convert to employment rate
dimnames(ER)[2] <- "ER"
chartSeries(ER,theme="white")
ER.data.frame <- as.data.frame(ER)
ER.data.frame$date <- ymd(rownames(ER.data.frame))
ER.time.series <- ts(ER.data.frame$ER, 
                     start = c(year(min(ER.data.frame$date)),month(min(ER.data.frame$date))),
                     end = c(year(max(ER.data.frame$date)),month(max(ER.data.frame$date))),
                     frequency=12)

test <- ER.time.series

#Test the function
cpi(test)
a <- cpi(test,1,5,2)
a[order(a$AIC),]