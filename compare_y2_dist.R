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
y = y2[y2>0]
library(survival)
m1=survreg(Surv(y)~1)
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
