library(dplyr)

setwd("C:\\Users\\dcries\\github\\bouts\\data")
pams <- read.csv('ISU_PAMS_2_minute_export.csv')

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


fcn <- function(data,skip,rejuv=10){
  n <- length(data)
  out <- matrix(NA,ncol=2,nrow=n)
  bout <- 1
  i <- 1
  while(i <= n){
    
    if((data[i] >= 3)){
      count <- 0
      nbout <- 1
      strt <- i
      i = i+1
      while((data[i] >= 3 | count <= skip) & (i <=n) & (count <= skip)){
        if(data[i] < 3){
          count = count+1
          i = i+1
        }
        else{
          i = i+1
        }
        
        if(i > rejuv){
          if((i/nbout >= rejuv+skip) & (all(data[(i-skip):(i-1)] > 3))){
            #print(count)
            #print(i)
            count <- 0
            nbout <- nbout+1
          }
        }

      }
      end <- i-1
      
      #eliminate possibility of last skip minutes being < 3 mets
      count2 <- 0
      while(data[end-count2] < 3){
        count2 <- count2+1
      }
      end2 <- end-count2
      
      
      i = i+1

      if((end2-strt>=9) & (sum(data[strt:end2])>=30)){
        out[bout,] <- c(strt,end2)
        bout <- bout+1
      }
    }
    
    else{
      i = i+1
    }
    
  }
  
  cout <- matrix(out[complete.cases(out),],ncol=2)
  
  m <- 0
  
  if(nrow(cout) > 0){
    for(i in 1:nrow(cout)){
      m <- m + sum(data[cout[i,1]:cout[i,2]]) - 30
    }
  }

  
  return(c(nrow(cout),m))
}