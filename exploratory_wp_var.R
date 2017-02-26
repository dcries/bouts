library(rjags)
library(ggplot2)
library(dplyr)

setwd("C:\\Users\\dcries\\github\\bouts\\data")
bouts <- read.csv("finalbouts.csv")

modelf <- "
model
{
  for(i in 1:n){
  y[i] ~ dpois(lambda[ind[i]] )
  }
  
  for(i in 1:n2){
    lambda[i] ~ dlnorm(mu, tau)
  }
  
  mu ~ dnorm(0,1/1000)
  tau ~ dgamma(1,1)

}
"
trues <- bouts %>% group_by(id) %>% filter(length(id)==2) %>% summarise(s=sd(nbouts))

meas2 <- unique(bouts$id)[which(table(bouts$id)==2)]
bouts2rep <- bouts[bouts$id %in% meas2,]
bouts2repc <- bouts2rep
#bouts2repc <- bouts2rep[complete.cases(bouts2rep),]
#bouts2repc <- bouts2repc %>% group_by(id) %>% filter(length(id)==2)
X <- bouts2repc[bouts2repc$rep==1,c("bmi","age","gender")]

datf <- list(y=bouts2repc$nbouts,
             ind=factor(bouts2repc$id,labels=1:length(unique(bouts2repc$id))),
             rep=bouts2repc$rep,
             n=nrow(bouts2repc),
             n2=length(unique(bouts2repc$id)),
             X=X,
             k=ncol(X))
mf = jags.model(textConnection(modelf), datf,n.adapt=1000,n.chains=3)
rf = coda.samples(mf, c("mu","tau","lambda"), n.iter=2000)

gelman.diag(rf[,c("mu","tau")])

#redefine trues so we are matching
trues <- bouts2repc %>% group_by(id) %>% filter(length(id)==2) %>% summarise(s=sd(nbouts),
                                                                             s2=sd(totalexcess))

lambda <- as.matrix(rf)[,3:(length(unique(bouts2repc$id))+2)]
nreps <- 1000
ind <- sample(1:nrow(lambda),nreps)
sims <- matrix(0,ncol=nrow(trues),nrow=nreps)
for(i in 1:nreps){
  for(j in 1:nrow(trues)){
    sims[i,j] <- sd(rpois(n=2,lambda = lambda[ind[i],j]))
  }
}

sims95 <- t(apply(sims,2,quantile,probs=c(0.025,0.5,0.975)))
full <- cbind(sims95,trues)
names(full)[1:3] <- c("low","med","high")
full$in_interval <- 0
full$in_interval[(full$s >= full$low) & (full$s <= full$high)] <- 1
full <- cbind(full,bouts2repc[bouts2repc$rep==1,c("bmi","age","gender")])

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

pvals <- rep(0,nrow(trues))
for(i in 1:nrow(trues)){
  pvals[i] <- sum(sims[,i] <= trues$s[i])/nreps
}

qplot(pvals) + geom_vline(xintercept=0.05,colour="red") + geom_vline(xintercept=0.95,colour="red") + theme_bw() + xlab("Posterior Predictive p-value")

sum(full$in_interval)/nrow(full) #0.87, undercovered ie. overdispersed for model
sum(full$in_interval[full$gender==1])/nrow(full[full$gender==1,])
sum(full$in_interval[full$gender==2])/nrow(full[full$gender==2,])
sum(full$in_interval[full$a==0])/nrow(full[full$a==0,])
sum(full$in_interval[full$a==1])/nrow(full[full$a==1,])
sum(full$in_interval[full$a==2])/nrow(full[full$a==2,])
sum(full$in_interval[full$a==3])/nrow(full[full$a==3,])
sum(full$in_interval[full$b==0])/nrow(full[full$b==0,])
sum(full$in_interval[full$b==1])/nrow(full[full$b==1,])
sum(full$in_interval[full$b==2])/nrow(full[full$b==2,])
sum(full$in_interval[full$b==3])/nrow(full[full$b==3,])

ggplot(data=full[1:100,]) + geom_pointrange(aes(x=as.factor(id),y=med,ymin=low,ymax=high,colour=as.factor(in_interval))) + 
  geom_point(aes(x=as.factor(id),y=s),colour="blue") + coord_flip()


qplot(data=full,y=s,x=(a),group=as.factor(a),geom="boxplot")
qplot(data=full,y=s,x=(b),group=as.factor(b),geom="boxplot")
qplot(data=full,y=s,x=gender,group=as.factor(gender),geom="boxplot")

qplot(data=full,y=s2,x=(a),group=as.factor(a),geom="boxplot")
qplot(data=full,y=s2,x=(b),group=as.factor(b),geom="boxplot")
qplot(data=full,y=s2,x=gender,group=as.factor(gender),geom="boxplot")

kruskal.test(s~a,data=full) #test variability depends on age
kruskal.test(s~b,data=full) #test variability depends on bmi
kruskal.test(s~gender,data=full) #test variability depends on gender

kruskal.test(s2~a,data=full) #test variability depends on age
kruskal.test(s2~b,data=full) #test variability depends on bmi
kruskal.test(s2~gender,data=full) #test variability depends on gender
