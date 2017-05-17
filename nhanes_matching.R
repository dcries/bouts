library(dplyr)
library(MatchIt)

setwd("C:\\Users\\dcries\\github\\bouts\\data")
bouts <- read.csv("finalbouts2rep.csv")

Za <- bouts %>% group_by(id) %>% filter(rep==1)
Za <- Za[,c("age","gender","bmi","smoke","education","black","hispanic")]
Za$education[Za$education <=3 ] <- 0
Za$education[Za$education >3 ] <- 1
Za$hispanic <- abs(Za$hispanic-2)

pams <- Za
pams$sample <- ('pams')
pams$id <- 1

nhanes <- read.csv("\\\\my.files.iastate.edu\\Users\\dcries\\Desktop\\research\\data\\NHANES_match.csv")
nhanes$sample <- "nhanes"

alldata <- rbind(pams,nhanes)
alldata$Group <- as.logical(alldata$sample == 'pams')

set.seed(1234)
match.it <- matchit(Group ~ age + gender + bmi+smoke+education+black+hispanic, data = alldata, method="nearest", ratio=1)
summary(match.it)

plot(match.it, type = 'jitter', interactive = FALSE)

df.match <- match.data(match.it)[1:ncol(alldata)]

sub <- df.match$id[df.match$Group==FALSE]

nhanes2 <- read.csv("\\\\my.files.iastate.edu\\Users\\dcries\\Desktop\\research\\data\\NHANES_final.csv")

nhanes3 <- nhanes2[nhanes2$id %in% sub,]

ids <- unique(nhanes3$id)
final <- matrix(0,ncol=ncol(nhanes3)+1,nrow=length(ids)*2)
for(i in 1:length(ids)){
  temp <- subset(nhanes3,id==ids[i])
  n <- nrow(temp)
  ind <- sample(1:n,2,FALSE)
  final[(2*i-1):(2*i),1:ncol(nhanes3)] <- as.matrix(temp[ind,])
  final[(2*i-1):(2*i),ncol(nhanes3)+1] <- 1:2
}

final <- data.frame(final)
names(final) <- c(names(nhanes3),"rep")

write.csv(final,file="\\\\my.files.iastate.edu\\Users\\dcries\\Desktop\\research\\data\\NHANES_use.csv")
write.csv(final,file="NHANES_use.csv")
