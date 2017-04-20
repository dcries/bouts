library(dplyr)
library(reshape)
library(ggplot2)
library(gridExtra)
library(MASS)

Rcpp::sourceCpp('C:/Users/dcries/github/bouts/bout_mcmc.cpp')
source('C:/Users/dcries/github/bouts/rgenpois.R')
source('C:/Users/dcries/github/bouts/pp_assess.R')
source('~/Documents/github/bouts/rgenpois.R')

setwd("C:\\Users\\dcries\\github\\bouts\\data")
setwd("~/Documents/github/bouts/data/")
bouts <- read.csv("finalbouts2rep.csv")

#Za <- bouts %>% group_by(id) %>% filter(rep==1) %>% select(age,gender,bmi,smoke,education,black,hispanic)
Za <- bouts %>% group_by(id) %>% filter(rep==1)
Za <- Za[,c("age","gender","bmi","smoke","education","black","hispanic")]
Za$education[Za$education <=3 ] <- 0
Za$education[Za$education >3 ] <- 1
Za$hispanic <- abs(Za$hispanic-2)

Za <- model.matrix(~age+as.factor(gender)+bmi+as.factor(smoke)+(education)+(black)+as.factor(hispanic),data=Za)
#Za[,2] <- scale(Za[,2])
#Za[,4] <- scale(Za[,4])

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
ncomp=3
data = list(Za=Za,Zb=Za,y1=y1,y2=y2)
init = list(currentbetay=rep(0,ncol(data$Za)+1),currentbetax=c(6.4,-.006,0.57,-.04,-.19,0,0,0,0.48),currentalpha=rep(0,ncol(data$Zb)),
            currentgamma=c(2.01,-0.012,0.578,-0.018,-0.011,0,0,0),currentsigma2y=0.95,currentsigma2x=6.73,
            currenteta=1.23,currentx1=rowMeans(data$y1)+0.1,currentx2=rowMeans(data$y2)+1,
            currentu=rep(0,nrow(data$y1)),gammatune=rep(0.00001,ncol(Za)),
            propa=1,propb=0.5,propx2=1/0.05,vx2=rep(10,nrow(Za)),
            x1propa=x1propa,x1propb=x1propb,
            currenttheta=rep(1,1),betaxtune=rep(0.0001,ncol(Za)+1), 
            propax2=1,propbx2=0.5,
            currentlambda=0.5,propl1=1,propl2=1,
            currentdelta=1,propd1=1,propd2=1,alphatune=rep(0.0000001,ncol(data$Zb)),
            currentb=matrix(0,nrow=nrow(data$y1),ncol=2),btune=c(0.001,0.001),
            currentSigmab=diag(2)*0.00000001,currentzeta=sample(1:ncomp,nrow(data$y1),TRUE),
            currentpi=rep(1/ncomp,ncomp),currentm=apply(y1,1,max)+2,
            currentsigma2b=0.0001,currentb2=matrix(0,nrow=nrow(data$y1),ncol=3),
            btune2=rep(0.001,3),currentSigmab2=diag(3)*0.01,currentb3=rep(0,nrow(data$y1)))

prior = list(mu0y2=rep(0,ncol(data$Za)+1),mu0x1=rep(0,ncol(Za)),mu0x2=rep(0,ncol(Za)+1),
             mu0a=rep(0,ncol(data$Zb)),V0y2=100*diag(ncol(data$Za)+1),V0x1=100*diag(ncol(Za)),
             V0x2=100*diag(ncol(Za)+1),V0a=100*diag(ncol(data$Zb)),a0eta=1,b0eta=1,
             a0x=1,b0x=1,a0y=1,b0y=1,
             a0theta=1,b0theta=1,
             a0l=1,b0l=1,
             a0delta=1,b0delta=1, d0=4, D0=diag(2),adirich=rep(1,ncomp),
             D02=diag(3))

mcmc = mcmc_2part_nci1(data=data,init=init,priors=prior,nrep=6000,burn=2000)
mcmc2 = mcmc_2part_nci2(data=data,init=init,priors=prior,nrep=6000,burn=2000)
mcmc3 = mcmc_2part_nci3(data=data,init=init,priors=prior,nrep=6000,burn=2000)
mcmc4 = mcmc_2part_nci4(data=data,init=init,priors=prior,nrep=6000,burn=2000)

mcmc2b = mcmc_2part_nci2b(data=data,init=init,priors=prior,nrep=6000,burn=2000)
mcmc2c = mcmc_2part_nci2c(data=data,init=init,priors=prior,nrep=6000,burn=2000)
mcmc2d = mcmc_2part_nci2d(data=data,init=init,priors=prior,nrep=6000,burn=2000)


#acceptance rates
apply(mcmc$gamma,2,function(x){return(length(unique(x))/length(x))}) 
#apply(mcmc$betax,2,function(x){return(length(unique(x))/length(x))}) 
#apply(mcmc$alpha,2,function(x){return(length(unique(x))/length(x))}) 
plot(apply(mcmc$latentx1,2,function(x){return((length(unique(x))-1)/length(x))}))
#plot(apply(mcmc$latentx2,2,function(x){return((length(unique(x))-1)/length(x))}))
plot(apply(mcmc$b1,2,function(x){return((length(unique(x))-1)/length(x))}))
(length(unique(mcmc$eta))-1)/length(mcmc$eta)
#(length(unique(mcmc$theta))-1)/length(mcmc$theta)
(length(unique(mcmc$lambda))-1)/length(mcmc$lambda)
#(length(unique(mcmc$delta))-1)/length(mcmc$delta)


hist(colMeans(mcmc$latentx1))
#hist(colMeans(mcmc$latentx2))

plot(mcmc$latentx1[,which.max(rowMeans(data$y1))])
#plot(mcmc$latentx2[,which.max(rowMeans(data$y2))])

plot(mcmc$gamma[,1],type="l")
plot(mcmc$gamma[,2],type="l")
plot(mcmc$gamma[,3],type="l")
plot(mcmc$gamma[,4],type="l")
plot(mcmc$gamma[,5],type="l")
plot(mcmc$gamma[,6],type="l")
plot(mcmc$gamma[,7],type="l")
plot(mcmc$gamma[,8],type="l")

plot(mcmc$eta,type="l")
#plot(mcmc$theta,type="l")
plot(mcmc$lambda,type="l")
#plot(mcmc$delta,type="l")


plot(mcmc$betay[,1],type="l")
plot(mcmc$betay[,2],type="l")
plot(mcmc$betay[,3],type="l")
plot(mcmc$betay[,4],type="l")
plot(mcmc$betay[,5],type="l")
plot(mcmc$betay[,6],type="l")
plot(mcmc$betay[,7],type="l")
plot(mcmc$betay[,8],type="l")
plot(mcmc$betay[,9],type="l")


# plot(mcmc$betax[,1],type="l")
# plot(mcmc$betax[,2],type="l")
# plot(mcmc$betax[,3],type="l")
# plot(mcmc$betax[,4],type="l")
# plot(mcmc$betax[,5],type="l")
# plot(mcmc$betax[,6],type="l")
# plot(mcmc$betax[,7],type="l")
# plot(mcmc$betax[,8],type="l")
# plot(mcmc$betax[,9],type="l")



plot(mcmc$sigma2b[,1],type="l")
plot(mcmc$sigma2b[,2],type="l")
plot(mcmc$sigma2b[,3],type="l")
plot(mcmc$corrb[,1],type="l")
plot(mcmc$corrb[,2],type="l")
plot(mcmc$corrb[,3],type="l")

plot(mcmc$sigma2y,type="l")


ind=which(apply(mcmc$latentx1,2,function(x){return((length(unique(x))-1)/length(x))})<0.1)
ind=which(y1rowmean>12)
i=ind[1];plot(mcmc$latentx1[,i]);y1[i,];summary(mcmc$latentx1[,i])
i=1;plot(mcmc$latentx2[,i]);y2[i,];summary(mcmc$latentx2[,i])


i=160;plot(mcmc$latentx1[,i]);y1[i,];summary(mcmc$latentx1[,i])


#-------------------------------------------------------------

assessln <- pp_assess(mcmc,data$Zb,200,1)
assessln2 <- pp_assess(mcmc2,data$Zb,200,2,burn=2000)
assessln3 <- pp_assess(mcmc3,data$Zb,200,3)
assessln4 <- pp_assess(mcmc4,data$Zb,200,4)
assessln2b <- pp_assess(mcmc2b,data$Zb,200,"2b")
assessln2c <- pp_assess(mcmc2c,data$Zb,200,"2c")
assessln2d <- pp_assess(mcmc2d,data$Zb,200,"2d")

y1zeroboth <- sum(rowSums(y1)==0)
y1zeroeither <- sum(apply(y1,1,function(x){return(!0%in%x)}))
y1ones <- sum(y1==1)
y1twos <- sum(y1==2)
y1meanwpsd <- mean(apply(y1,1,sd))
y1wprange <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
y1overallrange <- max(y1)-min(y1)
y2zeroboth <- sum(rowSums(y2)==0)
y2zeroeither <- sum(apply(y2,1,function(x){return(!0%in%x)}))
y2greaterthan <- sum(y1[,1]*30+y2[,1]>450/7)
y1y2regcoef <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff$diff))[2]
y1y2cor <- cor(c(y1[y1>0]),c(y2[y2>0]))
y2median <- median(c(data$y2))
y2q15 <- quantile(c(data$y2),probs=c(0.20))
y2q25 <- quantile(c(data$y2),probs=c(0.25))
y2q35 <- quantile(c(data$y2),probs=c(0.35))
y2q90 <- quantile(c(data$y2),probs=c(0.9))
y2daydiff <- mean(y2[,1]-y2[,2])

assess=assessln2d$out
q1 <- qplot(x=assess$y1zeroboth) + geom_vline(xintercept=y1zeroboth,colour="red") + theme_bw()
q2 <- qplot(x=assess$y1zeroeither) + geom_vline(xintercept=y1zeroeither,colour="red") + theme_bw()
q2b <- qplot(x=assess$y1ones) + geom_vline(xintercept=y1ones,colour="red") + theme_bw()
q2c <- qplot(x=assess$y1twos) + geom_vline(xintercept=y1twos,colour="red") + theme_bw()
q3 <- qplot(x=assess$y1meanwpsd) + geom_vline(xintercept=y1meanwpsd,colour="red") + theme_bw()
q4 <- qplot(x=assess$y1wprange,geom="bar") + geom_vline(xintercept=y1wprange,colour="red") + theme_bw()
#q5 <- qplot(x=assess$y1overallrange,geom="bar") + geom_vline(xintercept=y1overallrange,colour="red") + theme_bw()
q6 <- qplot(x=assess$y2greaterthan) + geom_vline(xintercept=y2greaterthan,colour="red") + theme_bw()
#q7 <- qplot(x=assess$y1y2regcoef) + geom_vline(xintercept=y1y2regcoef,colour="red") + theme_bw()
q7b <- qplot(x=assess$y1y2cor) + geom_vline(xintercept=y1y2cor,colour="red") + theme_bw()
q8 <- qplot(x=assess$y2median) + geom_vline(xintercept=y2median,colour="red") + theme_bw()
q9a <- qplot(x=assess$y2q15) + geom_vline(xintercept=y2q15,colour="red") + theme_bw()
q9b <- qplot(x=assess$y2q25) + geom_vline(xintercept=y2q25,colour="red") + theme_bw()
q9c <- qplot(x=assess$y2q35) + geom_vline(xintercept=y2q35,colour="red") + theme_bw()
q10 <- qplot(x=assess$y2q90) + geom_vline(xintercept=y2q90,colour="red") + theme_bw()
q11 <- qplot(x=assess$y2daydiff) + geom_vline(xintercept=y2daydiff,colour="red") + theme_bw()

grid.arrange(q1,q2,q2b,q2c,q3,q4,q6,q7b,q8,q9a,q9b,q9c,q10,q11)

sum(assess$y1zeroboth > y1zeroboth)/nrow(assess)
sum(assess$y1zeroeither > y1zeroeither)/nrow(assess)
sum(assess$y1meanwpsd > y1meanwpsd)/nrow(assess)
sum(assess$y1wprange > y1wprange)/nrow(assess)
#sum(assess$y1overallrange > y1overallrange)/nrow(assess)
sum(assess$y2greaterthan > y2greaterthan)/nrow(assess)
sum(assess$y1y2regcoef > y1y2regcoef)/nrow(assess)
sum(assess$y2median > y2median)/nrow(assess)
sum(assess$y2q25 > y2q25)/nrow(assess)
sum(assess$y2q90 > y2q90)/nrow(assess)
