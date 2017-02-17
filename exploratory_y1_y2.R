library(ggplot2)
library(dplyr)
library(reshape)
library(rjags)
library(scales)

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
# ggplot(data=y1,aes(x=nbouts)) + geom_bar() + theme_bw()
# ggplot(data=bybout,aes(x=totalp)) + geom_histogram(aes(y=..density..)) + theme_bw()
# ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + theme_bw()
# #by gender
# ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~gender) + theme_bw()
# ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~gender) + theme_bw()
# #by education
# ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~education)  + theme_bw()
# ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~education) + theme_bw()
# #smoke
# ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~smoke) + theme_bw()
# ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~smoke) + theme_bw()
# #employed
# ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~employed) + theme_bw()
# ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~employed) + theme_bw()
# #married
# ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~married) + theme_bw()
# ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~married) + theme_bw()
# #black
# ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~black) + theme_bw()
# ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~black) + theme_bw()
# #hispanic
# ggplot(data=y1,aes(x=nbouts)) + geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + facet_wrap(~hispanic) + theme_bw()
# ggplot(data=y2,aes(x=total)) + geom_histogram(aes(y=..density..)) + facet_wrap(~hispanic) + theme_bw()
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
bouts2rep$nbouts[bouts2rep$nbouts >= 10] <- 10 #cap nbouts at 20

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
    y[i] ~ dpois(lambda[ind[i]] + gamma[ind[i]] + delta[rep[i]])
  }

  for(i in 1:n2){
    lambda[i] ~ dgamma(alpha, beta)
    gamma[i] ~ dnorm(0,tau2g)
  }
  for(i in 1:2){
    delta[i] ~ dnorm(0,tau2d)
  }
  
  beta ~ dgamma(1,1)
  alpha ~ dgamma(1,1)
  tau2g ~ dgamma(1,1)
  tau2d ~ dgamma(1,1)
  sigmag <- 1/sqrt(tau2g)
  sigmad <- 1/sqrt(tau2d)

}
"

datf <- list(y=bouts$nbouts,
             ind=factor(bouts$id,labels=1:length(unique(bouts$id))),
             rep=bouts$rep,
             n=nrow(bouts),
             n2=length(unique(bouts$id)))
mf = jags.model(textConnection(modelf), datf,n.adapt=1000,n.chains=3)
rf = coda.samples(mf, c("alpha", "beta","sigmag", "sigmad" ,"lambda"), n.iter=2000)

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
kruskal.test(avgtotalexcess~nbouts,data=subset(bouts,nbouts>0&nbouts<11))


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


#contingency table for nbouts1 and nbouts2
table(trial1$bout,trial2$bout)
