library(dplyr)
library(reshape)
#library(ggplot2)
#library(gridExtra)

Rcpp::sourceCpp('/home/dcries/bouts/bout_mcmc_nci7.cpp')
source('/home/dcries/bouts/pp_assess.R')
source('/home/dcries/bouts/rgenpois.R')

#setwd("C:\\Users\\dcries\\github\\bouts\\data")
bouts <- read.csv("/home/dcries/bouts/data/finalbouts2rep.csv")
weights <- bouts %>% group_by(id) %>% filter(rep==1)
weights <- unlist(weights[,"B1BaseWeight"])
Za <- bouts %>% group_by(id) %>% filter(rep==1) %>% select(age,gender,bmi,smoke,education,black,hispanic)
Za$education[Za$education <=3 ] <- 0
Za$education[Za$education >3 ] <- 1
Za$hispanic <- abs(Za$hispanic-2)

Za <- model.matrix(~age+as.factor(gender)+bmi+as.factor(smoke)+(education)+(black)+as.factor(hispanic),data=Za)


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

Z1 <- rbind(data$Za,data$Za)
y1a <- c(y1)
m7 <- glm(y1a ~ as.matrix(Z1)+0,family=poisson)
valsg <- confint.default(m7,level=0.999)

Z = data.frame(rbind(data$Za[y2[,1]>0,],data$Za[y2[,2]>0,]))
names(Z) <- c("int","age","gender","bmi","smoke","education","black","hispanic")
y = y2[y2>0]
y1p = y1[y1>0]
m6 <- glm(y~as.matrix(Z)+log(y1p)+0,family=Gamma(link=power(lambda=1/2)))
valsb <- confint.default(m6,level=0.999)

data = list(Za=Za,Zb=Za,y1=y1,y2=y2)
init = list(currentbetay=valsb[,2],
            currentgamma=valsg[,2],currentsigma2y=0.95,currentsigma2x=6.73,
            currenteta=1.23,currentx1=rowMeans(data$y1)+0.1,currentx2=rowMeans(data$y2)+1,
            gammatune=rep(0.00000001,ncol(Za)),propa=1,propb=0.5,propx2=1/0.05,vx2=rep(10,nrow(Za)),
            x1propa=x1propa,x1propb=x1propb,betaxtune=c(1,rep(0.01,ncol(Za)-1),1), 
            propax2=1,propbx2=0.5,currentlambda=.1,propl1=1,propl2=1,
            propd1=1,propd2=1,currentb=matrix(0,nrow=nrow(data$y1),ncol=2),btune=c(0.001,0.001),
            currentSigmab=diag(2)*.01, currentsigma2b=0.01,currentphi=.39)

prior = list(mu0y2=rep(0,ncol(data$Za)+1),mu0x1=rep(0,ncol(Za)),mu0x2=rep(0,ncol(Za)+1),
             mu0a=rep(0,ncol(data$Zb)),V0y2=100*diag(ncol(data$Za)+1),V0x1=100*diag(ncol(Za)),
             V0x2=100*diag(ncol(Za)+1),V0a=100*diag(ncol(data$Zb)),a0eta=1,b0eta=1,
             a0x=1,b0x=1,a0y=1,b0y=1,
             a0theta=1,b0theta=1,
             a0l=1,b0l=1,
             a0delta=1,b0delta=1, d0=4, D0=diag(2))

mcmc = mcmc_2part_nci7(data=data,init=init,priors=prior,nrep=300000,burn=50000,thin=10)
save(mcmc,file="boutsmcmc3.RData")

assessln <- pp_assess(mcmc,data$Zb,1000,"7",y1,y2,weights,burn=0)

save(assessln,file="assessmcmc3.RData")