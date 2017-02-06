library(dplyr)

pams <- read.csv("~/Downloads/ISU_PAMS_2_minute_export.csv")

#
summary(pams$METS_per_minute_5.2[-1]-pams$METS_per_minute_5.2[-nrow(pams)])

plot(ecdf(pams$METS_per_minute_5.2[-1]-pams$METS_per_minute_5.2[-nrow(pams)]),xlim=c(-1,1))

diff <- pams %>% group_by(fileID) %>% summarise(nobs=length(METS_per_minute_5.2),
                                        q1 = quantile(METS_per_minute_5.2[-1]-METS_per_minute_5.2[-nobs],probs=0.01,na.rm=T),
                                        q5 = quantile(METS_per_minute_5.2[-1]-METS_per_minute_5.2[-nobs],probs=0.05,na.rm=T),
                                        q10 = quantile(METS_per_minute_5.2[-1]-METS_per_minute_5.2[-nobs],probs=0.1,na.rm=T),
                                        q90 = quantile(METS_per_minute_5.2[-1]-METS_per_minute_5.2[-nobs],probs=0.9,na.rm=T),
                                        q95 = quantile(METS_per_minute_5.2[-1]-METS_per_minute_5.2[-nobs],probs=0.95,na.rm=T),
                                        q99 = quantile(METS_per_minute_5.2[-1]-METS_per_minute_5.2[-nobs],probs=0.99,na.rm=T))


fcn <- function(data,skip){
  n <- length(data)
  out <- matrix(0,ncol=2,nrow=(n))
  bout <- 1
  i <- 1
  while(i < n){
    
    if(data[i] >= 3){
      count <- 0
      strt <- i
      i = i+1
      while(data[i] >= 3 | count < skip){
        if(data[i] < 3){
          count = count+1
          i = i+1
        }
        else{
          i = i+1
        }
      }
      end <- i-1
      i = i+1

      if((end-strt>=9) & (sum(data[strt:end])>=30)){
        out[bout,] < c(strt,end)
        bout <- bout+1
      }
    }
    
    else{
      i = i+1
    }
    
  }
  return(out)
}