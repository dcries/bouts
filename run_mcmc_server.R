library(dplyr)
library(reshape)
#library(ggplot2)
#library(gridExtra)

Rcpp::sourceCpp('/home/dcries/bouts/bout_mcmc.cpp')
#source('C:/Users/dcries/github/bouts/pp_assess.R')

#setwd("C:\\Users\\dcries\\github\\bouts\\data")
bouts <- read.csv("/home/dcries/bouts/data/finalbouts2rep.csv")

Za <- bouts %>% group_by(id) %>% filter(rep==1) %>% select(age,gender,bmi,smoke,education,black,hispanic)
Za$education[Za$education <=3 ] <- 0
Za$education[Za$education >3 ] <- 1
Za$hispanic <- abs(Za$hispanic-2)

Za <- model.matrix(~age+as.factor(gender)+bmi+as.factor(smoke)+(education)+(black)+as.factor(hispanic),data=Za)
Za[,2] <- scale(Za[,2])
Za[,4] <- scale(Za[,4])

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
init = list(currentbetay=c(0,0,1),currentbetax=c(6.4,-.006,0.57,-.04,-.19,0,0,0,0.48),currentalpha=rep(0,ncol(data$Zb)),
            currentgamma=c(2.01,-0.012,0.578,-0.018,-0.011,0,0,0),currentsigma2y=0.95,currentsigma2x=6.73,
            currenteta=1.23,currentx1=rowMeans(data$y1)+0.1,currentx2=rowMeans(data$y2)+1,
            currentu=rep(0,nrow(data$y1)),gammatune=rep(0.00001,ncol(Za)),
            propa=1,propb=0.5,propx2=1/0.05,vx2=rep(10,nrow(Za)),
            x1propa=x1propa,x1propb=x1propb,
            currenttheta=1,betaxtune=rep(0.0001,ncol(Za)+1), 
            propax2=1,propbx2=0.5,
            currentlambda=0.5,propl1=1,propl2=1,
            currentdelta=1,propd1=1,propd2=1,alphatune=rep(0.0000001,ncol(data$Zb)),
            currentb=matrix(0,nrow=nrow(data$y1),ncol=2),btune=c(0.001,0.001),
            currentSigmab=diag(2)*0.0001)

prior = list(mu0y2=rep(0,3),mu0x1=rep(0,ncol(Za)),mu0x2=rep(0,ncol(Za)+1),
             mu0a=rep(0,ncol(data$Zb)),V0y2=100*diag(3),V0x1=100*diag(ncol(Za)),
             V0x2=100*diag(ncol(Za)+1),V0a=100*diag(ncol(data$Zb)),a0eta=1,b0eta=1,
             a0x=1,b0x=1,a0y=1,b0y=1,
             a0theta=1,b0theta=1,
             a0l=1,b0l=1,
             a0delta=1,b0delta=1, d0=3, D0=diag(2))

mcmc = mcmc_2part_1(data=data,init=init,priors=prior,nrep=50000,burn=10000)
save(mcmc,file="boutsmcmc.RData")