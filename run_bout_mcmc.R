library(dplyr)
library(reshape)
library(MASS)
library(Hmisc)

Rcpp::sourceCpp('out_mcmc.cpp')
source("modelassess.R")
source('rgenpois.R')

bouts <- read.csv("finalbouts_impoajob.csv")
weights <- (bouts %>% group_by(id) %>% filter(rep==1))$rakedweights


Za <- bouts %>% group_by(id) %>% filter(rep==1) %>% dplyr::select(age,gender,smoke,education,black,hispanic,oajob)
Za <- model.matrix(~age+as.factor(gender)+as.factor(smoke)+(education)+(black)+as.factor(hispanic)+oajob,data=Za)

#put y1 and y2 in correct format for mcmc function
boutsmelt <- melt(bouts[,c("id","rep","nbouts","totalexcess")],id.vars=c("id","rep"))
boutscast <- cast(newm2,id~rep+variable)
y1=as.matrix(boutscast[,c(2,4)]);y2=as.matrix(boutscast[,c(3,5)])

#-----------------
#glm models to get starting values for beta and gamma
Z1 <- rbind(Za,Za)
y1a <- c(y1)
m1 <- glm(y1a ~ as.matrix(Z1)+0,family=poisson)
valsg <- confint.default(m1,level=0.9999)
Z = data.frame(rbind(Za[y2[,1]>0,],Za[y2[,2]>0,]))
names(Z) <- c("int","age","gender","smoke","education","black","hispanic","oajob")
y = y2[y2>0]
m2 <- glm(y~as.matrix(Z)+0,family=gaussian(link=log))
valsb <- confint.default(m2,level=0.9999)
#--------------

data = list(Za=Za,y1=y1,y2=y2)
#make y2 average excess MET-minutes
data$y2[data$y2[,1]>0,1] <- data$y2[data$y2[,1]>0,1]/data$y1[data$y2[,1]>0,1]
data$y2[data$y2[,2]>0,2] <- data$y2[data$y2[,2]>0,2]/data$y1[data$y2[,2]>0,2]


init = list(currentbetay=coef(m2),
            currentgamma=coef(m1),currentsigma2y=0.95,
            gammatune=rep(0.00000001,ncol(Za)),
            currentlambda=.5,proplambda1=1,proplambda2=1,
            currentb=matrix(0,nrow=nrow(data$y1),ncol=2),btune=c(0.0001,0.0001),
            currentSigmab=diag(2)*1, currentsigma2b=1)

prior = list(m_beta=rep(0,ncol(data$Za)),m_gamma=rep(0,ncol(Za)),
             V_beta=100*diag(ncol(data$Za)),V_gamma=100*diag(ncol(Za)),
             a_sigma2y=1,b_sigma2y=1,
             a0theta=1,b0theta=1,
             a_lambda=1,b_lambda=1,
             d0=4, D0=diag(2))

mcmc = mcmc_2part(data=data,init=init,priors=prior,nrep=3000,burn=2000,thin=1)

assessln <- pp_assess(mcmc,data$Zb,10,data$y1,data$y2,weights,burn=1)
