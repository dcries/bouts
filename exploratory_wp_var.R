library(rjags)
library(ggplot2)
library(dplyr)

setwd("C:\\Users\\dcries\\github\\bouts\\data")
bouts <- read.csv("finalbouts2rep.csv")

modelg <- "
model
{
  for(i in 1:n){
    #y[i] ~ dpois(lambda[ind[i]] )
    y[i] ~ dpois(lambda[i] )

    lambda[i] ~ dgamma(eta,theta[ind[i]])

    #y[i] ~ dpois(lambda2[ind[i]] )
  }
  
  for(i in 1:n2){
    mu[i] <- exp(inprod(X[i,],beta))#+epsilon[i]))
    theta[i] <- eta/mu[i]
    #lambda[i] ~ dgamma(eta,theta[i])
    #epsilon[i] ~ dnorm(0,taue)
    #lambda2[i] <- lambda[i]*(1-u[i])
    #u[i] ~ dbern(p[i])
    #logit(p[i]) <- inprod(X[i,],eta)
  }
  for(i in 1:k){
    beta[i] ~ dnorm(0,1/1000)
    #eta[i] ~ dnorm(0,1/1000)
  }
  eta ~ dgamma(1,1)
  #taue ~ dgamma(1,1)
  #sigmae <- 1/sqrt(taue)
  #p ~ dbeta(1,1)
}
"

modelln <- "
model
{
  for(i in 1:n){
    y[i] ~ dpois(lambda[ind[i]] )
  }
  
  for(i in 1:n2){
    lambda[i] ~ dlnorm(mu[i], tau)
    mu[i] <- inprod(X[i,],beta)
  }
  for(i in 1:k){
    beta[i] ~ dnorm(0,1/1000)
  }
  tau ~ dgamma(1,1)
  sigmau <- 1/sqrt(tau)
}
"
trues <- bouts %>% group_by(id) %>% summarise(s=(nbouts[2]-nbouts[1]),s2=sd(totalexcess),mean=mean(nbouts),total=nbouts[2]+nbouts[1])
X <- model.matrix(~bmi+age+factor(gender,labels=c(0,1),levels=c(1,2))+factor(smoke,labels=c(0,1),levels=c(1,2)),data=subset(bouts,rep==1))
#1 male, 0 female
#1 smoker, 0 non-smoker

datf <- list(y=bouts$nbouts,
             ind=factor(bouts$id,labels=1:length(unique(bouts$id))),
             rep=bouts$rep,
             n=nrow(bouts),
             n2=length(unique(bouts$id)),
             X=X,
             k=ncol(X))
mg = jags.model(textConnection(modelg), datf,n.adapt=1000,n.chains=3)
rg = coda.samples(mg, c("beta","eta","lambda"), n.iter=2000)
gelman.diag(rg[,c(paste0("beta[",1:ncol(X),"]"),"eta")])

#mln = jags.model(textConnection(modelln), datf,n.adapt=1000,n.chains=3)
#rln = coda.samples(mln, c("beta","sigmau","tau","lambda"), n.iter=2000)
#gelman.diag(rln[,c(paste0("beta[",1:ncol(X),"]"),"sigmau")])

#redefine trues so we are matching
#trues <- bouts2repc %>% group_by(id) %>% filter(length(id)==2) %>% summarise(s=(nbouts[2]-nbouts[1]),
#                                                                             s2=sd(totalexcess))
# #simulate from conditional dist poisson lognormal
# lambda <- as.matrix(rln)[,c(paste0("lambda[",1:nrow(trues),"]"))]
# nreps <- 1000
# ind <- sample(1:nrow(lambda),nreps)
# sims <- matrix(0,ncol=nrow(trues),nrow=nreps)
# for(i in 1:nreps){
#   for(j in 1:nrow(trues)){
#     sims[i,j] <- sd(rpois(n=2,lambda = lambda[ind[i],j]))
#   }
# }
# 
# #simulate from conditional dist poisson gamma
# lambdag <- as.matrix(rg)[,c(paste0("lambda[",1:nrow(trues),"]"))]
# nreps <- 1000
# ind <- sample(1:nrow(lambdag),nreps)
# simsg <- matrix(0,ncol=nrow(trues),nrow=nreps)
# for(i in 1:nreps){
#   for(j in 1:nrow(trues)){
#     simsg[i,j] <- sd(rpois(n=2,lambda = lambdag[ind[i],j]))
#   }
# }

#simulate from marginal for poisson lognormal
# n <- nrow(trues)
# beta <- as.matrix(rln)[,c(paste0("beta[",1:ncol(X),"]"))]
# sigmau <- unlist(rln[,"sigmau"])
# temp <- rep(0,2*n)
# m_med <- rep(0,nreps)
# m_range <- rep(0,nreps)
# m_0s <- rep(0,nreps)
# #m_max <- rep(0,nreps)
# m_count <- matrix(0,ncol=10,nrow=nreps)
# 
# within_diff <- matrix(0,nrow=nrow(trues),ncol=nreps)
# 
# for(i in 1:nreps){
#     lam <- rlnorm(n,X%*%(beta[ind[i],]),sigmau[ind[i]])
#   temp <- rpois(2*n,rep(lam,each=2))
#   tempmat <- matrix(temp,ncol=2,byrow=TRUE)
#   
#   within_diff[,i] <- apply(tempmat,1,function(x){return(x[2]-x[1])})
#   
#   m_med[i] <- median(temp)
#   m_range[i] <- range(temp)[2] - range(temp)[1]
#   m_0s[i] <- sum(temp==0)
#   for(j in 1:ncol(m_count)){
#     m_count[i,j] <- sum(temp==j)
#   }
# }


#simulate from marginal for poisson gamma
nreps <- 2000
n <- nrow(trues)
ind <- sample(1:nrow(lambda),nreps)
betag <- as.matrix(rg)[,c(paste0("beta[",1:ncol(X),"]"))]
eta <- unlist(rg[,"eta"])
tempg <- rep(0,2*n)
m_medg <- rep(0,nreps)
m_rangeg <- rep(0,nreps)
m_0sg <- rep(0,nreps)
#m_max <- rep(0,nreps)
m_countg <- matrix(0,ncol=10,nrow=nreps)
within_diffg <- matrix(0,nrow=nrow(trues),ncol=nreps)
within_mean <- matrix(0,nrow=nrow(trues),ncol=nreps)
within_sum <- matrix(0,nrow=nrow(trues),ncol=nreps)

#p <- unlist(rg[,"p"])
p <- as.matrix(rg)[,c(paste0("p[",1:nrow(X),"]"))]
sigmae <- unlist(rg[,"sigmae"])

for(i in 1:nreps){
  #lam <- rgamma(n,alpha[ind[i]],rate=beta[ind[i]])
  #lam <- rgamma(n,eta[ind[i]],rate=eta[ind[i]]/exp(X%*%(betag[ind[i],])))#+rnorm(n,0,sigmae[ind[i]])))
  #tempg <- rpois(2*n,rep(lam,each=2))#*rbinom(2*n,1,rep(1-p[ind[i],],each=2))
  lam <- rgamma(2*n,rep(eta[ind[i]],2*n),rate=rep(eta[ind[i]]/exp(X%*%(betag[ind[i],])),each=2))#+rnorm(n,0,sigmae[ind[i]])))
  tempg <- rpois(2*n,lam)#*rbinom(2*n,1,rep(1-p[ind[i],],each=2))
  
  tempmatg <- matrix(tempg,ncol=2,byrow=TRUE)
  
  within_diffg[,i] <- apply(tempmatg,1,function(x){return(x[2]-x[1])})
  within_mean[,i] <- rowMeans(tempmatg)
  within_sum[,i] <- rowSums(tempmatg)
  
  m_medg[i] <- median(tempg)
  m_rangeg[i] <- range(tempg)[2] - range(tempg)[1]
  m_0sg[i] <- sum(tempg==0)
  for(j in 1:ncol(m_count)){
    m_countg[i,j] <- sum(tempg==j)
  }
}

#----------------------------------------
# sims95 <- t(apply(sims,2,quantile,probs=c(0.025,0.5,0.975)))
# full <- cbind(sims95,trues)
# names(full)[1:3] <- c("low","med","high")
# full$in_interval <- 0
# full$in_interval[(full$s >= full$low) & (full$s <= full$high)] <- 1
# full <- cbind(full,bouts2repc[bouts2repc$rep==1,c("bmi","age","gender")])

full <- cbind(sims95,trues)
agegroup <- quantile(full$age,probs=c(0.25,0.5,0.75),na.rm=TRUE)
full$a <- 0
full$a[full$age > agegroup[1]] <- 1
full$a[full$age > agegroup[2]] <- 2
full$a[full$age > agegroup[3]] <- 3
bmigroup <- quantile(full$bmi,probs=c(0.25,0.5,0.75),na.rm=TRUE)
full$b <- 0
full$b[full$bmi > bmigroup[1]] <- 1
full$b[full$bmi > bmigroup[2]] <- 2
full$b[full$bmi > bmigroup[3]] <- 3

qplot(data=full,y=s,x=(a),group=as.factor(a),geom="boxplot")
qplot(data=full,y=s,x=(b),group=as.factor(b),geom="boxplot")
qplot(data=full,y=s,x=gender,group=as.factor(gender),geom="boxplot")

qplot(data=full,y=s2,x=(a),group=as.factor(a),geom="boxplot")
qplot(data=full,y=s2,x=(b),group=as.factor(b),geom="boxplot")
qplot(data=full,y=s2,x=gender,group=as.factor(gender),geom="boxplot")

# kruskal.test(s~a,data=full) #test variability depends on age
# kruskal.test(s~b,data=full) #test variability depends on bmi
# kruskal.test(s~gender,data=full) #test variability depends on gender
# 
# kruskal.test(s2~a,data=full) #test variability depends on age
# kruskal.test(s2~b,data=full) #test variability depends on bmi
# kruskal.test(s2~gender,data=full) #test variability depends on gender


#-------------------------------------------------------------
#check marginal for lognormal
# qplot(x=m_med) + geom_vline(xintercept=median(datf$y),colour="red")
# qplot(x=m_0s) + geom_vline(xintercept=sum(datf$y==0),colour="red")
# qplot(x=m_range) + geom_vline(xintercept=range(datf$y)[2]-range(datf$y)[1],colour="red")


#checking conditional for lognormal
# obs <- trues$s
# mean_obs <- mean(obs); sd_obs <- sd(obs); range_obs <- range(obs)[2]-range(obs)[1]
# mean_pp <- rowMeans(within_diff); sd_pp <- apply(within_diff,1,sd); range_pp <- apply(within_diff,1,function(x){return(range(x)[2]-range(x)[1])})
# 
# qplot(x=mean_pp) + geom_vline(xintercept=mean_obs,colour="red")
# qplot(x=sd_pp) + geom_vline(xintercept=sd_obs,colour="red")
# qplot(x=range_pp) + geom_vline(xintercept=range_obs,colour="red")

#check marginal for gamma
p1 <- qplot(x=m_medg) + geom_vline(xintercept=median(datf$y),colour="red") +xlab("median") + theme_bw()
p2 <- qplot(x=m_0sg) + geom_vline(xintercept=sum(datf$y==0),colour="red") +xlab("Number of Zeros")+ theme_bw()
p3 <- qplot(x=m_rangeg) + geom_vline(xintercept=range(datf$y)[2]-range(datf$y)[1],colour="red") +xlab("Range") + theme_bw()
grid.arrange(p1,p2,p3,nrow=1)

#checking conditional for gamma
obs <- trues$s
mean_obs <- mean(obs); sd_obs <- sd(obs); range_obs <- range(obs)[2]-range(obs)[1]
mean_ppg <- colMeans(within_diffg); sd_ppg <- apply(within_diffg,2,sd); range_ppg <- apply(within_diffg,2,function(x){return(range(x)[2]-range(x)[1])})

c1 <- qplot(x=mean_ppg) + geom_vline(xintercept=mean_obs,colour="red")+xlab("Mean") + theme_bw()
c2 <- qplot(x=sd_ppg) + geom_vline(xintercept=sd_obs,colour="red")+xlab("Standard Deviation") + theme_bw()
c3 <- qplot(x=range_ppg) + geom_vline(xintercept=range_obs,colour="red")+xlab("Range") + theme_bw()
grid.arrange(c1,c2,c3,nrow=1)


#posterior for beta
pg <- cbind(eta,betag)
mpg <- melt(pg)
ggplot(data=mpg,aes(x=value,y=..density..)) + geom_histogram() + facet_wrap(~X2,scales="free") + theme_bw()
