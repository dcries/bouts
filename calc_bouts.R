#this script reads the 8 minute by minute files and extracts only the "bout"
#data and eliminates the rest. It then pairs this data with the corresponding id
#so we can match demographic variables to these bout level data from the aggregate
#data. The function calc_bouts takes a vector of MET-mins by minute, how many 
#minutes an individual can be under 3 METs and we still continue counting the bout
# and the rejuvination time, ie. after x amount of minutes, we reset the skip to 0
#for example, an individual is >3 MEts for 8 min, 2 min below, then 8 more min >3
#METs, then 1 below, then 4 min >3, then 5 below. We would want to count it as a
#23 minute bout

#IMPORTANT!!!!!!
#a row in the outputted file 'bouts' that has 0 mets and 0 bout means that 
#individual never had 10+ min >3mets for that day, this is kind of just a placeholder


#the second part just combines the 8 output files into 1

library(lubridate)
#setwd("C:\\Users\\dcries\\github\\bouts\\data")
setwd('\\\\my.files.iastate.edu\\Users\\dcries\\Desktop\\research\\pams\\data')
#pams <- read.csv('ISU_PAMS_2_minute_export.csv')
finalpams <- read.csv('FinalPAMSDataSetNew.csv')

calc_bouts <- function(data,skip,rejuv=10){
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
  diff <- cout[,2]-cout[,1] + 1
  csdiff <- c(0,cumsum(diff))
  
  # m <- 0
  # 
  outmat <- matrix(NA,ncol=2,nrow=n)
  if(nrow(cout) > 0){
    for(i in 1:nrow(cout)){
      #m <- m + sum(data[cout[i,1]:cout[i,2]]) - 30
      outmat[(csdiff[i]+1):csdiff[i+1],1] <- data[cout[i,1]:cout[i,2]]
      outmat[(csdiff[i]+1):csdiff[i+1],2] <- i
    }
  }
  else{
    outmat[1,1] = outmat[1,2] = 0
  }

  
  return(outmat[complete.cases(outmat),])
}

set_day <- function(date){
  #must by POSIX class
  n <- length(date)
  dayn <- rep(0,n)
  i <- 1
  count <- 1
  while(i <= n){
    strtdate <- date[i]
    enddate <- date[i] + days(1)
    nmins <- sum((date >= strtdate) & (date < enddate))
    # if(nmins < 720){
    #   print("less than 1/2 day \n")
    #   if(nmins < 60){
    #     print("less than 60 minutes \n")
    #   }
    # }
    dayn[i:(i+nmins-1)] <- count
    
    check <- which(date >= enddate)
    i <- ifelse(length(check)>0, min(check), 9999999)
    count = count + 1
  }
  return(dayn)
}

nzeros <- rep(0,8)
bouts1 <- matrix(NA,nrow=nrow(pams),ncol=3)
for(i in 4:8){
  bouts1 <- matrix(NA,nrow=nrow(pams),ncol=3)
  
  pams <- read.csv(paste0('ISU_PAMS_',i,'_minute_export.csv'))
  fileid <- read.csv(paste0('ISUfileIDs_Q',i,'.csv'),header=FALSE)
  
  #rand#remove NAs?
  pams <- pams[complete.cases(pams),]
  
  pams$date <- ymd_hm(paste0(pams$Year,"-",pams$Month,"-",pams$Day, " ", pams$Hour, ":", pams$Minute))
  pams <- pams %>% group_by(fileID) %>% mutate(nday=set_day(date))
  #eliminate data for which <90% of day is recorded?
  cutoff <- 0.9*1440
  mins <- pams %>% group_by(fileID) %>% summarise(total1=sum(nday==1),
                                            total2=sum(nday==2),
                                            total3=sum(nday==3),
                                            total4=sum(nday==4))

  rm1 <- mins$fileID[mins$total1 < cutoff]
  rm2 <- mins$fileID[mins$total2 < cutoff]
  rm3 <- mins$fileID[mins$total3 < cutoff]
  rm4 <- mins$fileID[mins$total4 < cutoff]
  
  a=pams[ !(((pams$fileID %in% rm1)&(pams$nday==1)) | ((pams$fileID %in% rm2)&(pams$nday==2)) | ((pams$fileID %in% rm3)&(pams$nday==3)) | ((pams$fileID %in% rm4)&(pams$nday==4))) ,]
  
  for(j in 1:length(unique(pams$fileID))){
    temp <- matrix(calc_bouts(pams$METS_per_minute_5.2[pams$fileID==j],2,10),ncol=2)
    temp2 <- cbind(temp,rep(fileid[fileid[,1]==j,2],nrow(temp)))
    
    bouts1 <- rbind(bouts1,temp2)
    #print(j) not too slow
  }
  
  #remove NAs to get correct dimension
  bouts1 <- as.data.frame(bouts1[complete.cases(bouts1),])
  names(bouts1) <- c("mets","bout","id")
  
  for(j in 1:nrow(bouts1)){
    bouts1$age[j] <- finalpams$Age[bouts1$id[j]==finalpams$id][1]
    bouts1$gender[j] <- finalpams$Gender[bouts1$id[j]==finalpams$id][1]
    bouts1$bmi[j] <- finalpams$BMI[bouts1$id[j]==finalpams$id][1]
    bouts1$education[j] <- finalpams$Education[bouts1$id[j]==finalpams$id][1]
    bouts1$employed[j] <- finalpams$Employed[bouts1$id[j]==finalpams$id][1]
    bouts1$income[j] <- finalpams$Income[bouts1$id[j]==finalpams$id][1]
    bouts1$black[j] <- finalpams$Black_Ethnicity[bouts1$id[j]==finalpams$id][1]
    bouts1$hispanic[j] <- finalpams$Hispanic[bouts1$id[j]==finalpams$id][1]
    bouts1$smoke[j] <- finalpams$Smoker[bouts1$id[j]==finalpams$id][1]
    bouts1$occupation[j] <- finalpams$Occupation[bouts1$id[j]==finalpams$id][1]
    bouts1$married[j] <- finalpams$Married[bouts1$id[j]==finalpams$id][1]
    bouts1$weekend[j] <- finalpams$Weekend[bouts1$id[j]==finalpams$id][1]
    
    if(j%%1000==0){print(j)} #takes some time
  }
  
  write.csv(bouts1,file=paste0('bouts',i,'.csv'),row.names=FALSE)

}


#----------------------------------------------------------
#----------------------------------------------------------

bouts1 <- read.csv("bouts1.csv")
bouts2 <- read.csv("bouts2.csv")
bouts3 <- read.csv("bouts3.csv")
bouts4 <- read.csv("bouts4.csv")
bouts5 <- read.csv("bouts5.csv")
bouts6 <- read.csv("bouts6.csv")
bouts7 <- read.csv("bouts7.csv")
bouts8 <- read.csv("bouts8.csv")

bouts <- rbind(bouts1,bouts2,bouts4,bouts5,bouts6,bouts7,bouts8)
write.csv(bouts,file="bouts.csv",row.names=FALSE)

#----------------------------------------------------------