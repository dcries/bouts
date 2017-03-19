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

data = list(Za=Za,Zb=Za,y1=as.matrix(newc[,c(2,4)]),y2=as.matrix(newc[,c(3,5)]))
init = list(currentbetay=rep(0,3),currentbetax=rep(0,ncol(Za)+1),currentalpha=rep(0,3),
            currentgamma=rep(0,ncol(Za)),currentsigma2y=100,currentsigma2x=100,
            currenteta=2,currentx1=rowMeans(data$y1)+0.1,currentx2=rowMeans(data$y2)+1,
            currentu=rep(0,nrow(data$y1)),gammatune=rep(0.01,ncol(Za)),
            propa=0.4,propb=1,propx2=1/0.05,vx2=rep(10,nrow(Za)))

prior = list(mu0y2=rep(0,3),mu0x1=rep(0,ncol(Za)),mu0x2=rep(0,ncol(Za)+1),
             mu0a=rep(0,3),V0y2=100*diag(3),V0x1=100*diag(ncol(Za)),
             V0x2=100*diag(ncol(Za)+1),V0a=100*diag(3),a0eta=1,b0eta=1,
             a0x=1,b0x=1,a0y=1,b0y=1)

mcmc = mcmc_2part_1(data=data,init=init,priors=prior,nrep=2000,burn=1000)

#acceptance rates
apply(mcmc$gamma,2,function(x){return(length(unique(x))/length(x))})
plot(apply(mcmc$latentx1,2,function(x){return((length(unique(x))-1)/length(x))}))
plot(apply(mcmc$latentx2,2,function(x){return((length(unique(x))-1)/length(x))}))
length(unique(mcmc$eta))/length(mcmc$eta)

hist(colMeans(mcmc$latentx1))
hist(colMeans(mcmc$latentx2))

