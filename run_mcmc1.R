library(dplyr)
library(reshape)

Rcpp::sourceCpp('C:/Users/dcries/github/bouts/bout_mcmc.cpp')

setwd("C:\\Users\\dcries\\github\\bouts\\data")
bouts <- read.csv("finalbouts2rep.csv")

Za <- bouts %>% group_by(id) %>% filter(rep==1) %>% select(age,gender,bmi,smoke)
Za <- model.matrix(~age+as.factor(gender)+bmi+as.factor(smoke),data=Za)

new <- bouts[,c("id","rep","nbouts","totalexcess")]
newm <- melt(new,id.vars=c("id","rep"))
newc <- cast(newm,id~rep+variable)
x1tune=rowMeans(data$y1)
x1tune[x1tune==0] <- 1
y1rowmean=rowMeans(data$y1)
y1rowmean[y1rowmean==0] <- 1
y1rowvar=apply(data$y1,1,var)
y1rowvar[y1rowvar==0] <- 1
x1propb <- y1rowmean/y1rowvar
x1propa <- y1rowmean^2/y1rowvar
#mom estimator for et, alpha=407, beta=329
data = list(Za=Za,Zb=Za,y1=as.matrix(newc[,c(2,4)]),y2=as.matrix(newc[,c(3,5)]))
init = list(currentbetay=c(4.01,.06,.002),currentbetax=c(6.4,-.006,0.57,-.04,-.19,0.48),currentalpha=c(-.32,.47,.001),
            currentgamma=c(2.01,-0.012,0.578,-0.018,-0.011),currentsigma2y=0.835,currentsigma2x=1.73,
            currenteta=1.23,currentx1=rowMeans(data$y1)+0.1,currentx2=rowMeans(data$y2)+1,
            currentu=rep(0,nrow(data$y1)),gammatune=c(0.00001,0.00001,0.00001,0.00001,0.00001),
            propa=407,propb=329,propx2=1/0.05,vx2=rep(10,nrow(Za)),
            x1tune=x1tune,x1propa=x1propa,x1propb=x1propb)

prior = list(mu0y2=rep(0,3),mu0x1=rep(0,ncol(Za)),mu0x2=rep(0,ncol(Za)+1),
             mu0a=rep(0,3),V0y2=100*diag(3),V0x1=100*diag(ncol(Za)),
             V0x2=100*diag(ncol(Za)+1),V0a=100*diag(3),a0eta=1,b0eta=1,
             a0x=1,b0x=1,a0y=1,b0y=1)

mcmc = mcmc_2part_1(data=data,init=init,priors=prior,nrep=20000,burn=10000)

#acceptance rates
apply(mcmc$gamma,2,function(x){return(length(unique(x))/length(x))}) 
plot(apply(mcmc$latentx1,2,function(x){return((length(unique(x))-1)/length(x))}))
plot(apply(mcmc$latentx2,2,function(x){return((length(unique(x))-1)/length(x))}))
length(unique(mcmc$eta))/length(mcmc$eta)

hist(colMeans(mcmc$latentx1))
hist(colMeans(mcmc$latentx2))

plot(mcmc$latentx1[,which.max(rowMeans(data$y1))])
plot(mcmc$latentx2[,which.max(rowMeans(data$y2))])

plot(mcmc$gamma[,1])
plot(mcmc$gamma[,2])
plot(mcmc$gamma[,3])
plot(mcmc$gamma[,4])
plot(mcmc$gamma[,5])

plot(mcmc$eta)

plot(mcmc$betay[,1])
plot(mcmc$betay[,2])
plot(mcmc$betay[,3])

plot(mcmc$betax[,1])
plot(mcmc$betax[,2])
plot(mcmc$betax[,3])
plot(mcmc$betax[,4])
plot(mcmc$betax[,5])
plot(mcmc$betax[,6])

plot(mcmc$alpha[,1])
plot(mcmc$alpha[,2])
plot(mcmc$alpha[,3])

plot(mcmc$sigma2x)
plot(mcmc$sigma2y)
