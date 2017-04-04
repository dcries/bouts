library(dplyr)
library(reshape)
library(ggplot2)
library(gridExtra)

Rcpp::sourceCpp('C:/Users/dcries/github/bouts/bout_mcmc.cpp')
source('C:/Users/dcries/github/bouts/rgenpois.R')
source('C:/Users/dcries/github/bouts/pp_assess.R')

setwd("C:\\Users\\dcries\\github\\bouts\\data")
bouts <- read.csv("finalbouts2rep.csv")

Za <- bouts %>% group_by(id) %>% filter(rep==1) %>% select(age,gender,bmi,smoke)
Za <- model.matrix(~age+as.factor(gender)+bmi+as.factor(smoke),data=Za)

new <- bouts[,c("id","rep","nbouts","totalexcess")]
newm <- melt(new,id.vars=c("id","rep"))
newc <- cast(newm,id~rep+variable)
y1=as.matrix(newc[,c(2,4)]);y2=as.matrix(newc[,c(3,5)])
#x1tune=rowMeans(y1)
#x1tune[x1tune==0] <- 1
y1rowmean=rowMeans(y1)
y1rowmean[y1rowmean==0] <- 1
y1rowvar=apply(y1,1,var) #+ 1
y1rowvar[y1rowvar==0] <- 0.51
x1propb <- y1rowmean/y1rowvar
x1propa <- y1rowmean^2/y1rowvar
#mom estimator for et, alpha=407, beta=329
# muy 4.01,.06,.002
#currentalpha=c(-.32,0.47,.001)
data = list(Za=Za,Zb=Za,y1=y1,y2=y2)
init = list(currentbetay=c(0,0,1),currentbetax=c(6.4,-.006,0.57,-.04,-.19,0.48),currentalpha=rep(0,ncol(data$Zb)),
            currentgamma=c(2.01,-0.012,0.578,-0.018,-0.011),currentsigma2y=0.95,currentsigma2x=6.73,
            currenteta=1.23,currentx1=rowMeans(data$y1)+0.1,currentx2=rowMeans(data$y2)+1,
            currentu=rep(0,nrow(data$y1)),gammatune=c(0.00001,0.00001,0.00001,0.00001,0.00001),
            propa=1,propb=0.5,propx2=1/0.05,vx2=rep(10,nrow(Za)),
            x1propa=x1propa,x1propb=x1propb,
            currenttheta=1,betaxtune=c(0.0001,0.00001,0.00001,0.00001,0.00001,0.0001), 
            propax2=1,propbx2=0.5,
            currentlambda=0.5,propl1=1,propl2=1,
            currentdelta=1,propd1=1,propd2=1,alphatune=rep(0.0000001,ncol(data$Zb)),
            currentb=matrix(0,nrow=nrow(data$y1),ncol=3),btune=c(0.001,0.001,0.001),
            currentSigmab=diag(3)*0.0001)

prior = list(mu0y2=rep(0,3),mu0x1=rep(0,ncol(Za)),mu0x2=rep(0,ncol(Za)+1),
             mu0a=rep(0,ncol(data$Zb)),V0y2=100*diag(3),V0x1=100*diag(ncol(Za)),
             V0x2=100*diag(ncol(Za)+1),V0a=100*diag(ncol(data$Zb)),a0eta=1,b0eta=1,
             a0x=1,b0x=1,a0y=1,b0y=1,
             a0theta=1,b0theta=1,
             a0l=1,b0l=1,
             a0delta=1,b0delta=1, d0=14, D0=5*diag(3))

mcmc = mcmc_2part_1(data=data,init=init,priors=prior,nrep=3000,burn=1000)

#acceptance rates
apply(mcmc$gamma,2,function(x){return(length(unique(x))/length(x))}) 
apply(mcmc$betax,2,function(x){return(length(unique(x))/length(x))}) 
apply(mcmc$alpha,2,function(x){return(length(unique(x))/length(x))}) 
plot(apply(mcmc$latentx1,2,function(x){return((length(unique(x))-1)/length(x))}))
plot(apply(mcmc$latentx2,2,function(x){return((length(unique(x))-1)/length(x))}))
plot(apply(mcmc$b1,2,function(x){return((length(unique(x))-1)/length(x))}))
(length(unique(mcmc$eta))-1)/length(mcmc$eta)
(length(unique(mcmc$theta))-1)/length(mcmc$theta)
(length(unique(mcmc$lambda))-1)/length(mcmc$lambda)
(length(unique(mcmc$delta))-1)/length(mcmc$delta)


hist(colMeans(mcmc$latentx1))
hist(colMeans(mcmc$latentx2))

plot(mcmc$latentx1[,which.max(rowMeans(data$y1))])
plot(mcmc$latentx2[,which.max(rowMeans(data$y2))])

plot(mcmc$gamma[,1],type="l")
plot(mcmc$gamma[,2],type="l")
plot(mcmc$gamma[,3],type="l")
plot(mcmc$gamma[,4],type="l")
plot(mcmc$gamma[,5],type="l")

plot(mcmc$eta,type="l")
plot(mcmc$theta,type="l")
plot(mcmc$lambda,type="l")
plot(mcmc$delta,type="l")


plot(mcmc$betay[,1],type="l")
plot(mcmc$betay[,2],type="l")
plot(mcmc$betay[,3],type="l")

plot(mcmc$betax[,1],type="l")
plot(mcmc$betax[,2],type="l")
plot(mcmc$betax[,3],type="l")
plot(mcmc$betax[,4],type="l")
plot(mcmc$betax[,5],type="l")
plot(mcmc$betax[,6],type="l")

plot(mcmc$alpha[,1],type="l")
plot(mcmc$alpha[,2],type="l")
plot(mcmc$alpha[,3],type="l")
plot(mcmc$alpha[,4],type="l")
plot(mcmc$alpha[,5],type="l")

plot(mcmc$sigmab[,1],type="l")
plot(mcmc$sigmab[,2],type="l")
plot(mcmc$sigmab[,3],type="l")

plot(mcmc$corrb[,1],type="l")
plot(mcmc$corrb[,2],type="l")
plot(mcmc$corrb[,3],type="l")

plot(mcmc$sigma2x)
plot(mcmc$sigma2y,type="l")


ind=which(apply(mcmc$latentx1,2,function(x){return((length(unique(x))-1)/length(x))})<0.1)
ind=which(y1rowmean>12)
i=ind[1];plot(mcmc$latentx1[,i]);y1[i,];summary(mcmc$latentx1[,i])
i=1;plot(mcmc$latentx2[,i]);y2[i,];summary(mcmc$latentx2[,i])


i=160;plot(mcmc$latentx1[,i]);y1[i,];summary(mcmc$latentx1[,i])


#-------------------------------------------------------------

weekenddiff <- bouts %>% group_by(id) %>% summarise(diff=Weekend[1]-Weekend[2])
assess <- pp_assess(mcmc,data$Zb,weekenddiff$diff,200)

y1zeroboth <- sum(rowSums(y1)==0)
y1zeroeither <- sum(apply(y1,1,function(x){return(!0%in%x)}))
y1meanwpsd <- mean(apply(y1,1,sd))
y1wprange <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
y1overallrange <- max(y1)-min(y1)
y2zeroboth <- sum(rowSums(y2)==0)
y2zeroeither <- sum(apply(y2,1,function(x){return(!0%in%x)}))
y2greaterthan <- sum(y2>450/7)
y1y2regcoef <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff$diff))[2]

q1 <- qplot(x=assess$y1zeroboth) + geom_vline(xintercept=y1zeroboth,colour="red") + theme_bw()
q2 <- qplot(x=assess$y1zeroeither) + geom_vline(xintercept=y1zeroeither,colour="red") + theme_bw()
q3 <- qplot(x=assess$y1meanwpsd) + geom_vline(xintercept=y1meanwpsd,colour="red") + theme_bw()
q4 <- qplot(x=assess$y1wprange,geom="bar") + geom_vline(xintercept=y1wprange,colour="red") + theme_bw()
q5 <- qplot(x=assess$y1overallrange,geom="bar") + geom_vline(xintercept=y1overallrange,colour="red") + theme_bw()
q6 <- qplot(x=assess$y2greaterthan) + geom_vline(xintercept=y2greaterthan,colour="red") + theme_bw()
q7 <- qplot(x=assess$y1y2regcoef) + geom_vline(xintercept=y1y2regcoef,colour="red") + theme_bw()

grid.arrange(q1,q2,q3,q4,q5,q6,q7,nrow=3)

p1 <- sum(assess$y1zeroboth > y1zeroboth)/nrow(assess)
p2 <- sum(assess$y1zeroeither > y1zeroeither)/nrow(assess)
p3 <- sum(assess$y1meanwpsd > y1meanwpsd)/nrow(assess)
p4 <- sum(assess$y1wprange > y1wprange)/nrow(assess)
p5 <- sum(assess$y1overallrange > y1overallrange)/nrow(assess)
p6 <- sum(assess$y2greaterthan > y2greaterthan)/nrow(assess)
p7 <- sum(assess$y1y2regcoef > y1y2regcoef)/nrow(assess)
