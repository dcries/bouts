ll.weibull<-function(dat,par){
  a=(par[1])
  b=(par[2])
  
  ll=-sum(dweibull(dat,a,scale=b,log=T))
  
}

ll.gengamma<-function(dat,par){
  a1=(par[1])
  a2=(par[2])
  a3=par[3]
  ll=-sum(dgengamma.orig(dat,a1,a2,a3,log=T))
  
}
Z = data.frame(rbind(data$Za[y2[,1]>0,],data$Za[y2[,2]>0,]))
names(Z) <- c("int","age","gender","bmi","smoke","education","black","hispanic")
y = y2[y2>0]
y1p = y1[y1>0]
library(survival)
m1=survreg(Surv(y)~1)
m2=survreg(Surv(y)~age+gender+bmi+smoke+education+black+hispanic+log(y1p),data=Z)
m3=survreg(Surv(y)~log(y1p),data=Z)
m4=glm(y~age+gender+bmi+smoke+education+black+hispanic+log(y1p),data=Z,family=Gamma(link="log"))
m5=glm(y~age+gender+bmi+smoke+education+black+hispanic+log(y1p),data=Z,family=gaussian(link="log"))

a = 1/m1$scale
b = exp(coef(m1))
#shape 1 scale 2
par = c(a,b)
ks.test(y,"pweibull",a,b)
weibull.optim<-optim(par=par,fn=ll.weibull,dat=y)


alp <- mean(y)^2/var(y)
bet <- mean(y)/var(y)
ks.test(y,"pgamma",alp,bet)
ks.test(y,"plnorm",mean(log(y)),sd(log(y)))

gengamma.optim<-optim(par=c(3,2,2),fn=ll.weibull,dat=ym)


library(DEoptim)
ll.weibull<-function(par){
  a<-(par[1])
  b<-(par[2])
  
  ll=-sum(dweibull(y,a,scale=b,log=T))
  
}

lower <- c(0.001,0.001)
upper <- c(5,500)

outDEoptim <- DEoptim(ll.weibull, lower, upper,DEoptim.control(NP = 50,itermax = 200,trace=FALSE))
c(outDEoptim$optim$bestmem)


ll.gengamma<-function(par){
  a1=(par[1])
  a2=(par[2])
  a3=par[3]
  ll=-sum(dgengamma.orig(y,a1,a2,a3,log=T))
  
}

lower <- c(0.001,0.001,0.001)
upper <- c(500,500,500)

outDEoptim <- DEoptim(ll.gengamma, lower, upper,DEoptim.control(NP = 50,itermax = 200,trace=FALSE))
c(outDEoptim$optim$bestmem)
ks.test(y,"pgengamma.orig",0.55,53.83,2.29)

ggplot()+geom_histogram(aes(x=y,y=..density..))+geom_line(aes(x=1:2000,y=dweibull(1:2000,a,b)),colour="red")
ggplot()+geom_histogram(aes(x=y,y=..density..))+geom_line(aes(x=1:2000,y=dgengamma.orig(1:2000,0.55,53.83,2.29)),colour="red")

ysim = rgengamma.orig(10000,0.553,53.837,2.295)
ysim2 = rweibull(10000,par[1],par[2])
ysim3 = rgamma(10000,alp,bet)
ysim4 = rlnorm(10000,mean(log(y)),sd(log(y)))
plot(ecdf(y))
lines(ecdf(ysim),col="red")
lines(ecdf(ysim2),col="blue")
lines(ecdf(ysim3),col="green")
lines(ecdf(ysim4),col="brown")

df <- data.frame(gengamma=ysim,weibull=ysim2,gamma=ysim3,lnormal=ysim4)
mdf <- melt(df)
real <- data.frame(cbind(rep("real",length(y)),y))
names(real) <- names(mdf)
mdf <- rbind(mdf,real)
mdf$value <- as.numeric(mdf$value)

ggplot() + geom_density(data=mdf,aes(x=value,colour=variable)) + theme_bw() + xlim(c(0,2000)) + scale_colour_manual(values=c("red","blue","green","yellow","black"))


#-----------------------------------
#make power link function
library(survival)
library(rjags)
Z = data.frame(rbind(data$Za[y2[,1]>0,],data$Za[y2[,2]>0,]))
names(Z) <- c("int","age","gender","bmi","smoke","education","black","hispanic")
y = y2[y2>0]
y1p = y1[y1>0]

m2=survreg(Surv(y)~age+gender+bmi+smoke+education+black+hispanic+log(y1p),data=Z)
m3=survreg.old(Surv(y)~age+gender+bmi+smoke+education+black+hispanic+log(y1p),data=Z,link=power(lambda=1/2))

m4=glm(y~age+gender+bmi+smoke+education+black+hispanic+log(y1p),data=Z,family=Gamma(link="log"))

m6 <- glm(y~as.matrix(Z)+log(y1p)+0,family=Gamma(link=power(lambda=1/2)))



model <- "model{
for(i in 1:n){
  mu[i] <- (inprod(Z[i,],beta[]))^2
  #y[i] ~ dgamma(delta,delta/mu[i])
  y[i] ~ dweibull(delta,mu[i])
}
for(i in 1:k){
  beta[i] ~ dnorm(0.01,1/100)
}
delta ~ dgamma(1,1)
}
"
dattp <- list(y=y,n=length(y),Z=cbind(as.matrix(Z),log(y1p)),k=ncol(Z)+1)
mtp = jags.model(textConnection(model), dattp,n.adapt=1000,n.chains=3)
rtp = coda.samples(mtp, c("delta","beta"), n.iter=1000)
