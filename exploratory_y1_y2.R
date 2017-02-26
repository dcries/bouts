library(ggplot2)
library(dplyr)
library(reshape)
library(rjags)
library(scales)
library(rcompanion)
library(reshape)

setwd("C:\\Users\\dcries\\github\\bouts\\data")
boutsraw <- read.csv("bouts.csv")
bouts <- read.csv("finalbouts.csv")

bybout = boutsraw %>% group_by(id,rep,bout) %>% summarise(totalp = sum(mets),
                                                    totalpadj = sum(mets)-30,
                                                    age=age[1],
                                                    gender=gender[1],
                                                    bmi=bmi[1],
                                                    education=education[1],
                                                    employed=employed[1],
                                                    income=income[1],
                                                    black=black[1],
                                                    hispanic=hispanic[1],
                                                    smoke=smoke[1],
                                                    occupation=occupation[1],
                                                    married=married[1])
bybout$totalpadj[bybout$totalpadj==-30] <- 0
bybout$boutadj <- bybout$bout
bybout$boutadj[bybout$bout >= 10] <- 10
# y2 <- bybout %>% group_by(id,rep) %>% summarise(nbouts = max(bout) ,
#                                             total = sum(totalp),
#                                             totalexcess = sum(totalp2),
#                                              age=age[1],
#                                              gender=gender[1],
#                                              bmi=bmi[1],
#                                              education=education[1],
#                                              employed=employed[1],
#                                              income=income[1],
#                                              black=black[1],
#                                              hispanic=hispanic[1],
#                                              smoke=smoke[1],
#                                              occupation=occupation[1],
#                                              married=married[1])
#y2 <- y2[y2$total!=0,]
# a=y2 %>% group_by(id) %>% summarise(total2=sum(total))
# sum(a$total2==0)
# 
ggplot(data=bouts,aes(x=nbouts)) + geom_bar() + theme_bw()
ggplot(data=bybout,aes(x=totalp)) + geom_histogram(aes(y=..density..)) + theme_bw()
ggplot(data=bouts,aes(x=total)) + geom_histogram(aes(y=..density..)) + theme_bw()
#by gender
ggplot(data=bouts,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~gender) + theme_bw()
ggplot(data=bouts,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~gender) + theme_bw()
#by education
ggplot(data=bouts,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~education)  + theme_bw()
ggplot(data=bouts,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~education) + theme_bw()
#smoke
ggplot(data=bouts,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~smoke) + theme_bw()
ggplot(data=bouts,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~smoke) + theme_bw()
#employed
ggplot(data=bouts,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~employed) + theme_bw()
ggplot(data=bouts,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~employed) + theme_bw()
#married
ggplot(data=bouts,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~married) + theme_bw()
ggplot(data=bouts,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~married) + theme_bw()
#black
ggplot(data=bouts,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~black) + theme_bw()
ggplot(data=bouts,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~black) + theme_bw()
#hispanic
ggplot(data=bouts,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~hispanic) + theme_bw()
ggplot(data=bouts,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~hispanic) + theme_bw()
# 
# 
# #all
# plot(ecdf(bybout$totalp))
# plot(ecdf(y2$total))
# 
# #rep 1
# plot(ecdf(bybout$totalp[bybout$rep==1]))
# plot(ecdf(bybout$totalp[bybout$rep==2]))
# # diff
# aa <- bybout %>% group_by(id,rep) %>% summarise(s=sum(totalp))
# aaa <- aa %>% group_by(id) %>% filter(length((id))==2) %>% summarise(diff=s[2]-s[1])
# plot(ecdf(aaa$diff))
# 
# 
# #y1 vs y2
# plot(y1$nbouts,y2$total)
# plot(y1$nbouts,y1$age)
# plot(y1$nbouts,y1$bmi)
# plot(y2$total,y2$bmi)
# plot(y2$total,y2$age)


#--------------------------------------------------#
#--------------------------------------------------#

#bouts <- read.csv("finalbouts.csv")

meas2 <- unique(bouts$id)[which(table(bouts$id)==2)]
bouts2rep <- bouts[bouts$id %in% meas2,]
bouts2rep$nbouts[bouts2rep$nbouts >= 11] <- 11 #cap nbouts at 20

heat <- bouts2rep %>% group_by(id) %>% summarise(nbouts1=nbouts[1],
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
# p <- ggplot(heat,aes(nbouts1,nbouts2)) 
# p + stat_bin2d(aes(fill = ..count..), binwidth = c(1,1)) + theme_bw()
# p + stat_density(aes(fill = ..count..), geom = "raster")
# 
# breaks <- c(0, 2, 4,6,8,10,20,40,60)
# p + stat_bin2d(aes(fill = ..count..), binwidth = c(1,1))  +
#   scale_fill_gradientn(colours = rainbow(7),
#                        breaks = breaks, labels = format(breaks)) + theme_bw()
# 
# p + stat_bin2d(aes(fill = cut(..count..,c(0,2,5,10,50,100))), binwidth = c(1,1))
# p + stat_bin2d(aes(fill = ..count..), binwidth = c(1,1)) + scale_fill_gradient2(breaks=c(1,3,6,10,30,50))
#----------------------------------------------------#
bouts$nbouts2 <- bouts$nbouts
bouts$nbouts2[bouts$nbouts >= 5] <- 5

qplot(data=bouts,x=nbouts2,y=total,group=nbouts2,geom="boxplot")
qplot(data=bouts,x=nbouts2,y=totalexcess,group=nbouts2,geom="boxplot")

#--------------
#for poisson-gamma mixture, mle of alpha and beta
loglik <- function(par,y){
  n <- length(y)
  sumy <- sum(y)
  alpha <- par[1]
  beta <- n*1/((alpha+sumy)/alpha-1) #par[2]
  num <- alpha*log(beta) + lgamma(sumy+alpha)
  denom <- lgamma(alpha) + sum(lfactorial(y)) + (alpha+sumy)*log(beta+n)
  ll <- num-denom
  return(-ll)
}
pars <- optim(c(.9,0.3),loglik,y=bouts$nbouts)
pars <- optimize(loglik,c(0.0001,100000),y=bouts$nbouts)

#

lambda <- rgamma(10000,l^2,rate=l)
simy <- rpois(10000,lambda)
pg_l_1 <- nrow(bouts)*table(simy)/length(simy)
pg_l_1 <- c(pg_l_1,rep(0,max(bouts$nbouts)+1-length(pg_l_1)))

alpha2 <- l^2/5 #mean=l, var=1
beta2 <- l/5
lambda2 <- rgamma(10000,alpha2,rate=beta2)
simy2 <- rpois(10000,lambda2)
pg_l_5 <- nrow(bouts)*table(simy2)/length(simy2)
pg_l_5 <- c(pg_l_5,rep(0,max(bouts$nbouts)+1-length(pg_l_5)))

modelf <- "
model
{
  for(i in 1:n){
    y[i] ~ dpois(lambda[ind[i]] )#+ gamma[ind[i]] + delta[rep[i]])
  }

  for(i in 1:n2){
    lambda[i] ~ dgamma(alpha, beta)
    #gamma[i] ~ dnorm(0,tau2g)
  }
  for(i in 1:2){
    #delta[i] ~ dnorm(0,tau2d)
  }
  
  beta ~ dgamma(1,1)
  alpha ~ dgamma(1,1)
  #tau2g ~ dgamma(1,1)
  #tau2d ~ dgamma(1,1)
  #sigmag <- 1/sqrt(tau2g)
  #sigmad <- 1/sqrt(tau2d)

}
"

datf <- list(y=bouts$nbouts,
             ind=factor(bouts$id,labels=1:length(unique(bouts$id))),
             rep=bouts$rep,
             n=nrow(bouts),
             n2=length(unique(bouts$id)))
mf = jags.model(textConnection(modelf), datf,n.adapt=1000,n.chains=3)
rf = coda.samples(mf, c("alpha", "beta","lambda"), n.iter=2000)

alpha3 <- 0.984 #jags 1.150
beta3 <- 0.3033 # 0.3517
lambda3 <- rgamma(10000,alpha3,rate=beta3)
simy3 <- rpois(10000,lambda3)
pg_l_j <- nrow(bouts)*table(simy3)/length(simy3)
pg_l_j <- c(pg_l_j,rep(0,max(bouts$nbouts)+1-length(pg_l_j)))

#-------------------------

#fit of poisson model to nbouts
obs <- rep(0,max(bouts$nbouts)+1)
for(i in 0:max(bouts$nbouts)){
  obs[i+1] <- sum(bouts$nbouts==i)
}

l <- mean(bouts$nbouts)
lv <- nrow(bouts)*dpois(0:max(bouts$nbouts),l)

df <- data.frame(obs=obs,pois=lv,pg_l_j=pg_l_j,id=0:max(bouts$nbouts))
mdf <- melt(df,id.vars="id")

ggplot(data=mdf, aes(x=id,y=value,fill=variable)) + 
  geom_bar(stat="identity",alpha=0.5,position="dodge") + theme_bw() + xlim(c(-1,20))

#---------------------------
#Does min/bout increase with # of bouts?
qplot(data=bybout,x=boutadj,y=totalpadj,group=boutadj,geom="boxplot")

#--------------------------
#compare excess minutes time1 vs time2
 aa <- bybout %>% group_by(id,rep) %>% summarise(s=sum(totalp),
                                                 sadj=sum(totalpadj))
 aaa <- aa %>% group_by(id) %>% filter(length((id))==2) %>% summarise(diff=s[2]-s[1],
                                                                      diffadj=sadj[2]-sadj[1])

qplot(x=bouts2rep$totalexcess[bouts2rep$rep==1],y=bouts2rep$totalexcess[bouts2rep$rep==2]) + geom_abline(slope = 1,intercept = 0)
qplot(data=aaa,x=diff)
qplot(data=aaa,x=diffadj)


#------------------------------------------
#nbouts of weekday vs weekend
#qplot(data=bouts,x=nbouts,geom="bar",facets=~Weekend)
ggplot(data=subset(bouts,!is.na(Weekend)),aes(x=nbouts,y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + 
  geom_bar() + facet_wrap(~Weekend)

ks.test(bouts$nbouts[bouts$Weekend==0],bouts$nbouts[bouts$Weekend==1])
summary(lm(nbouts~Weekend,data=bouts)) #not good fit
summary(lm(nbouts~Weekend,data=subset(bouts,nbouts <= 20)))#not good fit
summary(glm(nbouts~Weekend,family="poisson",data=bouts))#not good fit
summary(glm(nbouts~Weekend,family=quasipoisson(link = "log"),data=subset(bouts,nbouts <= 20)))

summary(bouts$nbouts[bouts$Weekend==0])
summary(bouts$nbouts[bouts$Weekend==1])

View(bouts[which(is.na(bouts$Weekend)),])

#------------------------------
#excess minutes of weekday vs weekend
ggplot(data=bouts,aes(x=totalexcess,y=..density..)) + 
  geom_histogram() + facet_wrap(~Weekend)

ks.test(bouts$totalexcess[bouts$Weekend==0],bouts$totalexcess[bouts$Weekend==1])
summary(lm(totalexcess~Weekend,data=bouts)) #not good fit
(lm(log(totalexcess+.01)~Weekend,data=bouts))
m2 <- (glm((totalexcess+0.01)~Weekend,data=bouts,family="Gamma"))

summary(bouts$totalexcess[bouts$Weekend==0])
summary(bouts$totalexcess[bouts$Weekend==1])

#-----------------------------------------
#tests
#difference in totalexcess min from trial 1 to trial 2
trial1 <- bouts %>% group_by(id) %>% filter(length((id))==2, rep==1)
trial2 <- bouts %>% group_by(id) %>% filter(length((id))==2, rep==2)

t.test(trial1$totalexcess,trial2$totalexcess,paired=TRUE)
wilcox.test(trial1$totalexcess,trial2$totalexcess,paired=TRUE)
qqnorm(trial1$totalexcess-trial2$totalexcess);qqline(trial1$totalexcess-trial2$totalexcess)
shapiro.test(trial1$totalexcess-trial2$totalexcess)

#test distributional differences
ks.test(bouts$totalexcess[bouts$rep==1],bouts$totalexcess[bouts$rep==2])

#difference in nbouts min from trial 1 to trial 2
t.test(trial1$nbouts,trial2$nbouts,paired=TRUE)
wilcox.test(trial1$nbouts,trial2$nbouts,paired=TRUE)
qqnorm(trial1$nbouts-trial2$nbouts);qqline(trial1$nbouts-trial2$nbouts)
shapiro.test(trial1$nbouts-trial2$nbouts)

#test distributional differences
ks.test(bouts$nbouts[bouts$rep==1],bouts$nbouts[bouts$rep==2])

#nbouts of weekday vs weekend
ks.test(bouts$nbouts[bouts$Weekend==0],bouts$nbouts[bouts$Weekend==1])
#excess minutes of weekday vs weekend
ks.test(bouts$totalexcess[bouts$Weekend==0],bouts$totalexcess[bouts$Weekend==1])

#diff in weekday vs weekend for those with both
week1weekend1 <- bouts %>% group_by(id) %>% filter(length((id))==2, sum(Weekend)==1, rep==1)
week1weekend2 <- bouts %>% group_by(id) %>% filter(length((id))==2, sum(Weekend)==1, rep==2)

t.test(week1weekend1$totalexcess,week1weekend2$totalexcess,paired=TRUE)
wilcox.test(week1weekend1$totalexcess,week1weekend2$totalexcess,paired=TRUE)
qqnorm(week1weekend1$totalexcess-week1weekend2$totalexcess);qqline(week1weekend1$totalexcess-week1weekend2$totalexcess)
shapiro.test(week1weekend1$totalexcess-week1weekend2$totalexcess)

#
#test for difference in avg total excess mins by number of bouts, assumes normality which doesn't hold
m1 <- lm(avgtotalexcess~as.factor(nbouts),data=subset(bouts,nbouts>0))
anova(m1)
m1b <- lm(avgtotalexcess~as.factor(nbouts),data=subset(bouts,nbouts>0&nbouts<11))
anova(m1b)
m1lm <- lm(avgtotalexcess~(nbouts),data=subset(bouts,nbouts>0))
m1lmb <- lm(avgtotalexcess~(nbouts),data=subset(bouts,nbouts>0&nbouts<11))

#nonparametric test since normality doesn't hold
kruskal.test(avgtotalexcess~nbouts,data=subset(bouts,nbouts>0))
kruskal.test(avgtotalexcess~nbouts,data=subset(bouts,nbouts>0&nbouts<13))


#same test but taking log
qqnorm(log(bouts$avgtotalexcess[bouts$nbouts>0]));qqline(log(bouts$avgtotalexcess[bouts$nbouts>0]))
shapiro.test(log(bouts$avgtotalexcess[bouts$nbouts>0]))

m1l <- lm(log(avgtotalexcess)~as.factor(nbouts),data=subset(bouts,nbouts>0))
anova(m1l)
m1bl <- lm(log(avgtotalexcess)~as.factor(nbouts),data=subset(bouts,nbouts>0&nbouts<11))
anova(m1bl)
m1lm <- lm(log(avgtotalexcess)~(nbouts),data=subset(bouts,nbouts>0))
summary(m1lm)
m1lmb <- lm(log(avgtotalexcess)~(nbouts),data=subset(bouts,nbouts>0&nbouts<11))
summary(m1lmb)

#test for diff in total excess mins by bout number, assumes normality which doesn't hold
m2 <- lm(totalpadj~as.factor(bout),data=subset(bybout,bout>0))
anova(m2)
m2b <- lm(totalpadj~as.factor(bout),data=subset(bybout,bout>0&bout<11))
anova(m2b)
m2lm <- lm(totalpadj~(bout),data=subset(bybout,bout>0))
m2lmb <- lm(totalpadj~(bout),data=subset(bybout,bout>0&bout<11))

#nonparametric test since normality doesn't hold
kruskal.test(totalpadj~bout,data=subset(bybout,bout>0))
kruskal.test(totalpadj~bout,data=subset(bybout,bout>0&bout<11))


#contingency table for nbouts1 and nbouts2, but this is paired data
chisq.test(table(trial1$nbouts,trial2$nbouts))

#------------------------
#
ct0 <- matrix(table(trial1$nbouts,trial2$nbouts),nrow=21,byrow=T)
ct1 <- matrix(table(trial1$nbouts,trial2$nbouts)[1:11,1:11],nrow=11,byrow=T)
ct1diff <- data.frame(obs2=rep(0:9,10),obs1=rep(0:9,each=10),count=c(ct1-t(ct1)))


ggplot(data=ct1diff,aes(x=as.factor(obs1),y=as.factor(obs2))) + geom_tile(aes(fill=count)) + 
  geom_text(aes(label=count)) + scale_fill_gradientn (
    colours = colorRampPalette(c("blue" ,"white", "red"))(120),
    values = c(0,0.5,1))

ggplot(data=ct1diff,aes(x=as.factor(obs1),y=as.factor(obs2))) + geom_tile(aes(fill=abs(count))) + 
  geom_text(aes(label=abs(count))) + scale_fill_gradient(low = "white", high = "red") 
  

#-----------------------------
#McNemar test
mcnemar.test(ct1)
mcnemar.test(ct1[1:6,1:6])

# stuart maxwell test
library(irr)
stuart.maxwell.mh(ct1[1:6,1:6])



nominalSymmetryTest(ct1[1:6,1:6],MonteCarlo=TRUE,ntrial=100000,method = "fdr")


#--------------------------------------
#permutation test
#H_0 : sum abs lower.tri-upper.tri  = 0, ie. trial 1 and trial 2 exchangeable
# H_a : sum != 0, not exchangeable
obs <- sum(abs(ct1[lower.tri(ct1)]-ct1[upper.tri(ct1)]))
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
nsim <- 10000
permutei <- rep(0,nsim)
nboutmat <- cbind(trial1$nbouts,trial2$nbouts)
newmat <- matrix(0,nrow=nrow(nboutmat),ncol=ncol(nboutmat))
for(i in 1:nsim){
  for(j in 1:nrow(nboutmat)){
    s <- sample(1:2,2)
    newmat[j,s[1]] <- nboutmat[j,1]
    newmat[j,s[2]] <- nboutmat[j,2]
  }
  cti <- matrix(table(newmat[,1],newmat[,2])[1:10,1:10],ncol=10,byrow=T)
  #mat <- ct1[sample(nrow(ct1)),sample(ncol(ct1))]
  permutei[i] <- sum(abs(cti[lower.tri(cti)]-cti[upper.tri(cti)]))
}

qplot(x=permutei) + geom_vline(xintercept=obs,col="red")

#----------------------------------
#checking exchangeability for y2 total excess minutes
#or like this? instead permute obs from trial 1 and trial 2?
#obsy2 <- sum(trial1$totalexcess-trial2$totalexcess)#cor(trial1$totalexcess,trial2$totalexcess)
obsy2 <- coef(lm(trial1$totalexcess~trial2$totalexcess))[2]
nsim <- 10000
permutey2 <- rep(0,nsim)
nboutmaty2 <- cbind(trial1$totalexcess,trial2$totalexcess)
newmaty2 <- matrix(0,nrow=nrow(nboutmaty2),ncol=ncol(nboutmaty2))
for(i in 1:nsim){
  for(j in 1:nrow(nboutmaty2)){
    s <- sample(1:2,2)
    newmaty2[j,s[1]] <- nboutmaty2[j,1]
    newmaty2[j,s[2]] <- nboutmaty2[j,2]
  }
  
  #permutey2[i] <- sum(newmaty2[,1]-newmaty2[,2])#cor(newmaty2)[1,2]
  permutey2[i] <- coef(lm(newmaty2[,1]~newmaty2[,2]))[2]#cor(newmaty2)[1,2]
  
}

qplot(x=permutey2) + geom_vline(xintercept=obsy2,col="red") + theme_bw()

#----------------------------------------------------#
# long
#-----------------------------------------------------
#checking exchangeability for y2 total excess minutes conditional on nbouts
#nbouts is taken to be mean of two observations
permute_y2 <- function(data,nsim){
  data <- as.matrix(data)
  obs <- coef(lm(data[,1]~data[,2]))[2]
  permute <- rep(0,nsim)
  #nboutmat <- as.matrix(y2check[y2check$mnbouts==0.5,c("totalexcess1","totalexcess2")])
  newmat <- matrix(0,nrow=nrow(data),ncol=ncol(data))
  for(i in 1:nsim){
    for(j in 1:nrow(data)){
      s <- sample(1:2,2)
      newmat[j,s[1]] <- data[j,1]
      newmat[j,s[2]] <- data[j,2]
    }
    permute[i] <- coef(lm(newmat[,1]~newmat[,2]))[2]#cor(newmaty2)[1,2]
  }
  return(permute)
}

y2check <- bouts2rep %>% group_by(id) %>% summarise(mnbouts=mean(nbouts),
                                         totalexcess1=totalexcess[1],
                                         totalexcess2=totalexcess[2]
                                         )
nsim <- 1000

#mean nbouts = 0.5
obsy2_05 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==0.5)))[2]
permutey2_05 <- permute_y2(y2check[y2check$mnbouts==0.5,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_05) + geom_vline(xintercept=obsy2_05,col="red") + theme_bw()

#mean nbouts = 1
obsy2_1 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==1)))[2]
permutey2_1 <- permute_y2(y2check[y2check$mnbouts==1,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_1) + geom_vline(xintercept=obsy2_1,col="red") + theme_bw()

#mean nbouts = 1.5
obsy2_15 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==1.5)))[2]
permutey2_15 <- permute_y2(y2check[y2check$mnbouts==1.5,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_15) + geom_vline(xintercept=obsy2_15,col="red") + theme_bw()

#mean nbouts = 2
obsy2_2 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==2)))[2]
permutey2_2 <- permute_y2(y2check[y2check$mnbouts==2,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_2) + geom_vline(xintercept=obsy2_2,col="red") + theme_bw()

#mean nbouts = 2.5
obsy2_25 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==2.5)))[2]
permutey2_25 <- permute_y2(y2check[y2check$mnbouts==2.5,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_25) + geom_vline(xintercept=obsy2_25,col="red") + theme_bw()

#mean nbouts = 3
obsy2_3 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==3)))[2]
permutey2_3 <- permute_y2(y2check[y2check$mnbouts==3,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_3) + geom_vline(xintercept=obsy2_3,col="red") + theme_bw()

#mean nbouts = 3.5
obsy2_35 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==3.5)))[2]
permutey2_35 <- permute_y2(y2check[y2check$mnbouts==3.5,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_35) + geom_vline(xintercept=obsy2_35,col="red") + theme_bw()

#mean nbouts = 4
obsy2_4 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==4)))[2]
permutey2_4 <- permute_y2(y2check[y2check$mnbouts==4,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_4) + geom_vline(xintercept=obsy2_4,col="red") + theme_bw()

#mean nbouts = 4.5
obsy2_45 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==4.5)))[2]
permutey2_45 <- permute_y2(y2check[y2check$mnbouts==4.5,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_45) + geom_vline(xintercept=obsy2_45,col="red") + theme_bw()

#mean nbouts = 5
obsy2_5 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==5)))[2]
permutey2_5 <- permute_y2(y2check[y2check$mnbouts==5,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_5) + geom_vline(xintercept=obsy2_5,col="red") + theme_bw()

#mean nbouts = 5.5
obsy2_55 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==5.5)))[2]
permutey2_55 <- permute_y2(y2check[y2check$mnbouts==5.5,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_55) + geom_vline(xintercept=obsy2_55,col="red") + theme_bw()

#mean nbouts = 6
obsy2_6 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==6)))[2]
permutey2_6 <- permute_y2(y2check[y2check$mnbouts==6,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_6) + geom_vline(xintercept=obsy2_6,col="red") + theme_bw()

#mean nbouts = 6.5
obsy2_65 <- coef(lm(totalexcess1~totalexcess2,data=subset(y2check,mnbouts==6.5)))[2]
permutey2_65 <- permute_y2(y2check[y2check$mnbouts==6.5,c("totalexcess1","totalexcess2")],nsim)
#qplot(x=permutey2_65) + geom_vline(xintercept=obsy2_65,col="red") + theme_bw()


permutey2_df <- cbind(permutey2_05,permutey2_1,permutey2_15,permutey2_2,
                      permutey2_25,permutey2_3,permutey2_35,permutey2_4,
                      permutey2_45,permutey2_5,permutey2_55,permutey2_6,
                      permutey2_65)

truey2 <- c(obsy2_05,obsy2_1,obsy2_15,obsy2_2,obsy2_25,
            obsy2_3,obsy2_35,obsy2_4,obsy2_45,obsy2_5,
            obsy2_55,obsy2_6,obsy2_65)

permutey2_pvals <- rep(0,ncol(permutey2_df))
for(i in 1:ncol(permutey2_df)){
  permutey2_pvals[i] <- sum(permutey2_df[,i]<truey2[i])/nsim
}


permutey2_mdf <- melt(permutey2_df)
permutey2_mdf$true <- rep(truey2,each=nsim)

ggplot(data=permutey2_mdf,aes(x=value,group=X2)) + geom_histogram() + 
  geom_vline(aes(xintercept=true),colour="red") + facet_wrap(~X2,scales="free") + theme_bw()

#-----------------------------------------------------

#----------------------------------------------------------#
#check exchangeability for y2 (and weekend effect) using waynes idea
#
bouts2rep <- bouts[bouts$id %in% meas2,]
a <- melt(bouts2rep[,c("id","rep","nbouts","totalexcess","Weekend")],id.vars=c("id","rep"))
boutscast <- cast(data=a,id~rep+variable,value.var="value")

names(boutscast)[2:7] <- c("nbouts1","totalexcess1","weekend1","nbouts2","totalexcess2","weekend2")
boutscast$nboutsdiff <- boutscast$nbouts1-boutscast$nbouts2
boutscast$totalexcessdiff <- boutscast$totalexcess1-boutscast$totalexcess2
boutscast$weekenddiff <- boutscast$weekend1-boutscast$weekend2

y2exch <- lm(totalexcessdiff~nboutsdiff+weekenddiff,data=boutscast)
summary(y2exch)
y2exch2 <- lm(totalexcessdiff~nboutsdiff+weekenddiff,data=subset(boutscast,nbouts1<15&nbouts2<15))
summary(y2exch2)
y2exch3 <- lm(totalexcessdiff~nboutsdiff+weekenddiff,data=subset(boutscast,abs(nboutsdiff)<11))
summary(y2exch3)

library(stargazer)
stargazer(y2exch)

#------------------------------------------------
qplot(data=bouts,nbouts,geom="bar",facets=~rep)
dendf <- data.frame(nbouts1=trial1$nbouts,nbouts2=trial2$nbouts)
ggplot(data=dendf,aes(x=nbouts1,y=nbouts2)) + geom_point() + geom_density2d()
ggplot(data=dendf,aes(x=nbouts1,y=nbouts2)) + 
  stat_density2d(aes(fill=..density..),geom="tile",contour=FALSE) + 
  scale_fill_gradient(low="white",high="red") + xlim(c(0,10)) + ylim(c(0,10))
