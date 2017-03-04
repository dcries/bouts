library(ggplot2)
library(dplyr)
library(reshape)
library(rjags)
library(scales)
library(rcompanion)
library(reshape)


setwd("C:\\Users\\dcries\\github\\bouts\\data")
#boutsraw <- read.csv("bouts.csv")
#bouts <- read.csv("finalbouts.csv") #excludes > 2000 y2
#boutsall <- read.csv("finalboutsfull.csv") #everyone
bouts <- read.csv("finalbouts2rep.csv")
pams <- read.csv("FinalPAMSDataSetNew.csv")

#remove individuals with only 1 rep, and without basic demographics, bmi,age,gender,smoker
#bouts <- bouts %>% group_by(id) %>% filter(length(id)==2, !is.na(age), !is.na(bmi), !is.na(gender), !is.na(smoke))
#write.csv(bouts,file="finalbouts2rep.csv",row.names=FALSE)

#---------------------------------------------------------------
heat <- bouts %>% group_by(id) %>% summarise(nbouts1=nbouts[1],
                                                 nbouts2=nbouts[2],
                                                 out=paste0(nbouts[1]," ",nbouts[2]),
                                                 diff=nbouts[1]-nbouts[2])

a1 <- table(sort(heat$out))

count <- rep(0,121)
for(i in 0:10){
  for(j in 0:10){
    count[i*11+j+1] <- sum(heat$nbouts1==i & heat$nbouts2==j)
  }
}

heatdf <- data.frame(nbouts1=rep(0:10,each=11),nbouts2=rep(0:10,11),count=count)
qn <- quantile(heatdf$count,probs=c(0.01,0.1,0.25,0.5,0.75,0.9,0.99))
qn01 <- rescale(c(qn, range(heatdf$count))) 

ggplot(data=heatdf,aes(x=nbouts1,y=nbouts2)) + geom_tile(aes(fill=count)) + 
  geom_text(aes(label=count)) + #+ scale_fill_gradient(low = "white", high = "red") 
  scale_fill_gradientn (
    colours = colorRampPalette(c( "white", "darkred"))(120),
    values = c(0,qn01[2:7],1))
#------------------------------------------------------------------------

#-----------------------------------------
#tests
#difference in totalexcess min from trial 1 to trial 2
trial1 <- bouts %>% group_by(id) %>% filter(rep==1) %>% arrange(id)
trial2 <- bouts %>% group_by(id) %>% filter(rep==2) %>% arrange(id)

#difference in total excess min from trial 1 to trial 2
t.test(trial1$totalexcess,trial2$totalexcess,paired=TRUE)
wilcox.test(trial1$totalexcess,trial2$totalexcess,paired=TRUE)

#difference in nbouts min from trial 1 to trial 2
t.test(trial1$nbouts,trial2$nbouts,paired=TRUE)
wilcox.test(trial1$nbouts,trial2$nbouts,paired=TRUE)

#diff in weekday vs weekend for those with both
week1weekend1 <- bouts %>% group_by(id) %>% filter(sum(Weekend)==1, Weekend==0) %>% arrange(id)
week1weekend2 <- bouts %>% group_by(id) %>% filter(sum(Weekend)==1, Weekend==1) %>% arrange(id)

t.test(week1weekend1$totalexcess,week1weekend2$totalexcess,paired=TRUE)
wilcox.test(week1weekend1$totalexcess,week1weekend2$totalexcess,paired=TRUE)

t.test(week1weekend1$nbouts,week1weekend2$nbouts,paired=TRUE)
wilcox.test(week1weekend1$nbouts,week1weekend2$nbouts,paired=TRUE)

#contingency talbe
ct1 <- matrix(table(trial1$nbouts,trial2$nbouts)[1:11,1:11],nrow=11,byrow=T)
#McNemar test
mcnemar.test(ct1)
mcnemar.test(ct1[1:6,1:6])

#--------------------------------------
#permutation test
#H_0 : sum abs lower.tri-upper.tri  = 0, ie. trial 1 and trial 2 exchangeable
# H_a : sum != 0, not exchangeable
obs <- sum(abs(ct1[lower.tri(ct1)]-ct1[upper.tri(ct1)]))
obs <- t.test(trial1$nbouts,trial2$nbouts,paired=TRUE)$statistic
obs <- mean(trial1$nbouts-trial2$nbouts)
# 
# nsim <- 10000
# permute <- rep(0,nsim)
# vals <- as.numeric(ct1)
# for(i in 1:nsim){
#   mat <- ct1[sample(nrow(ct1)),sample(ncol(ct1))]
#   permute[i] <- sum(abs(mat[lower.tri(mat)]-mat[upper.tri(mat)]))
# }
# 
# qplot(x=permute) + geom_vline(xintercept=obs,col="red")

#or like this? instead permute obs from trial 1 and trial 2?
nsim <- 1000
permutei <- rep(0,nsim)
nboutmat <- cbind(trial1$nbouts,trial2$nbouts)
newmat <- matrix(0,nrow=nrow(nboutmat),ncol=ncol(nboutmat))
for(i in 1:nsim){
  for(j in 1:nrow(nboutmat)){
    s <- sample(1:2,2)
    newmat[j,s[1]] <- nboutmat[j,1]
    newmat[j,s[2]] <- nboutmat[j,2]
  }
  cti <- matrix(table(newmat[,1],newmat[,2])[1:nrow(ct1),1:nrow(ct1)],ncol=nrow(ct1),byrow=T)
  mat <- ct1[sample(nrow(ct1)),sample(ncol(ct1))]
  #permutei[i] <- sum(abs(cti[lower.tri(cti)]-cti[upper.tri(cti)]))
  #permutei[i] <- t.test(newmat[,1],newmat[,2],paired=TRUE)$statistic
  permutei[i] <- mean(newmat[,1]-newmat[,2])
  
}

qplot(x=permutei) + geom_vline(xintercept=obs,col="red")

#----------------------------------


#----------------------------------------------------------#
#check exchangeability for y2 (and weekend effect) using waynes idea
#

a <- melt(bouts[,c("id","rep","nbouts","totalexcess","Weekend")],id.vars=c("id","rep"))
boutscast <- cast(data=a,id~rep+variable,value.var="value")

names(boutscast)[2:7] <- c("nbouts1","totalexcess1","weekend1","nbouts2","totalexcess2","weekend2")
boutscast$nboutsdiff <- boutscast$nbouts1-boutscast$nbouts2
boutscast$totalexcessdiff <- boutscast$totalexcess1-boutscast$totalexcess2
boutscast$weekenddiff <- boutscast$weekend1-boutscast$weekend2

y2exch <- lm(totalexcessdiff~nboutsdiff+weekenddiff,data=boutscast)
summary(y2exch)
