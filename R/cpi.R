#'AEDForecasting:  
#' Package to incorporate change point analysis in ARIMA forecasting
#' @param myts a time series object
#' @param  startChangePoint  a positive integer for minimum number of changepoints
#' @param  endChangePoint  a positive integer for maximum number of change points. If 0 then  only startChangePoint number of change points will be entered. Should be either 0 or greater than startChangePoint and if so the algorithm will loop through all values inbetween subject to step
#' @param  step  an integer to step through loop of change points
#' @param  num Bump model number (see below)
#' @param  cpmeth changepoint method. Default is BinSeg. See cpa package for details
#' @param  CPpenalty default is SIC. See cpa package for details

#'@return The model number, number of change points, and the information criteria values
#'
#' The bump models are described below. 

cpi <- function(myts, startChangePoint = 1, endChangePoint = 0, step = 1, num=8, cpmeth='BinSeg', CPpenalty="SIC") {
  
  if (is.ts(myts) == FALSE) {
    message("Error: First parameter must be a time series")
  } else if (!((startChangePoint > 0) && (startChangePoint %% 1 == 0))) {
    message("Error: Starting change point must be a positive integer")
  } else if (!((endChangePoint > 0) && (endChangePoint %% 1 == 0)) && (endChangePoint != 0)){
    
    message("Error: Ending change point must be a positive integer")
  } else if ((endChangePoint <= startChangePoint) && (endChangePoint != 0)) {
    message("Error: Ending change point must be greater than starting change point")
  } else {
 
    DF2 <- 0
  
  
    cpiP2 <- function(myts, n=startChangePoint, cpmeth.=cpmeth, CPpenalty.=CPpenalty, num.=num) {
      
      DF <- data.frame('Model'=rep("", num), 'No.ChangePoints'=rep(n, num),
                       'AIC' = rep(0,num),stringsAsFactors = FALSE)
      for (k in 1:num) {
      
        m.bin <- suppressWarnings(changepoint::cpt.mean(myts,penalty=CPpenalty,method=cpmeth,Q=n))
        aveTS=mean(myts)
        temp1 <- c()
        temp2 <- c()
        cpiv <- matrix(0,ncol = n, nrow = length(myts))
        for (i in 1:n) {
          
          for (j in 1:m.bin@cpts[i]) {
            cpiv[j,i] <- 0 #for the ith interv var fills cpiv with 0s until the ith cpa row.
          }
          
          temp1 <- length(myts) - m.bin@cpts[i]
          
          for (j in 1:temp1) {
            temp2 <- j + m.bin@cpts[i]
            if(k == 1) {
              cpiv[temp2,i] <- 1 - (1-(log(temp2)/temp2)) 
              #DF[1,1] <- "1 - (1 - log(changepoint/changepoint)"
              DF[1,1] <- "model 1"
            }#model 1
            
            else if (k==2) {
              cpiv[temp2,i] <- aveTS*(1 - (temp2-m.bin@cpts[i])/(1 +(temp2-m.bin@cpts[i]))) 
              DF[2,1] <- "model 2"
            }#model 2
            
            else if (k==3) {
              cpiv[temp2,i] <- aveTS*(1-(1 /(temp2-m.bin@cpts[i])))
              DF[3,1] <- "model 3"
            } #model 3
            
            else if (k==4) {
              cpiv[temp2,i] <- aveTS*(1-(1 /(temp2-m.bin@cpts[i])))
              DF[4,1] <- "model 4"
            } #model 4
            
            else if (k==5) {
              cpiv[temp2,i] <- aveTS*exp(-((temp2-m.bin@cpts[i])^2)/((length(myts)^2)*2.5))
              DF[5,1] <- "model 5"
            } #model 5
            
            else if (k==6) {
              cpiv[temp2,i] <- aveTS*(1-(1 /(temp2-m.bin@cpts[i])))
              DF[6,1] <- "model 6"
            } #model 6
            
            else if (k==7) {
              cpiv[temp2,i] <- 1 - (1-(1/log(temp2)))
              DF[7,1] <- "model 7"
            } #model 7
            
            else if (k==8) {
              cpiv[temp2,i] <- 1
              DF[8,1] <- "model 8"
            } 
          }
        }
        xxx=suppressWarnings(forecast::auto.arima(myts,xreg = cpiv))
        #message=paste(n,k, sep="  ")
        #print(message)
        #print(xxx$aic)
        DF[k,3] <- xxx$aic
      }
    
    return(DF)
    }
    
    suppressWarnings(if (endChangePoint == 0) {
      DF2 <- cpiP2(myts)
    } else {
      for (z in seq(startChangePoint,endChangePoint,step)) {
        if (DF2 == 0) {
          DF2 <- cpiP2(myts,n=z)
        } else {
          DFtemp <- cpiP2(myts,n=z)
          DF2 <- merge(DF2, DFtemp, all = TRUE)
        }
      }
    })
    
    return(DF2)
  }
}
