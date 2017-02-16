library(dplyr)
setwd("C:\\Users\\dcries\\github\\bouts\\data")
bouts <- read.csv("finalbouts.csv")
FinalPAMSDataSetNew <- read.csv("FinalPAMSDataSetNew.csv")

bouts$id2 <- bouts$id*(bouts$rep+1.2345)#*(bouts$age+1.2345)
#bouts$num <- 1:nrow(bouts)
FinalPAMSDataSetNew$id2 <- (FinalPAMSDataSetNew$id)*(FinalPAMSDataSetNew$Trial+1.2345)#*(FinalPAMSDataSetNew$Age+1.2345)

boutsweekend= left_join(bouts,FinalPAMSDataSetNew[,c("id2","Weekend")],by="id2")

write.csv(boutsweekend,file="finalbouts.csv",row.names=FALSE)
