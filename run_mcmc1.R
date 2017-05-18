library(dplyr)
library(reshape)
library(ggplot2)
library(gridExtra)
library(MASS)
library(scales)

Rcpp::sourceCpp('C:/Users/dcries/github/bouts/bout_mcmc.cpp')
source('C:/Users/dcries/github/bouts/rgenpois.R')
source('C:/Users/dcries/github/bouts/pp_assess.R')
source('~/Documents/github/bouts/rgenpois.R')

setwd("C:\\Users\\dcries\\github\\bouts\\data")
setwd("~/Documents/github/bouts/data/")
bouts <- read.csv("finalbouts2rep.csv")
pams <- read.csv("FinalPAMSDataSetNew.csv")

#Za <- bouts %>% group_by(id) %>% filter(rep==1) %>% select(age,gender,bmi,smoke,education,black,hispanic)
weights <- bouts %>% group_by(id) %>% filter(rep==1)
weights <- unlist(weights[,"B1BaseWeight"])
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
init = list(currentbetay=c(1,0,0,0,0,0,0,0,0),currentbetax=c(6.4,-.006,0.57,-.04,-.19,0,0,0,0.48),currentalpha=rep(0,ncol(data$Zb)),
            currentgamma=c(2.01,-0.012,0.578,-0.018,-0.011,0,0,0),currentsigma2y=0.95,currentsigma2x=6.73,
            currenteta=1.23,currentx1=rowMeans(data$y1)+0.1,currentx2=rowMeans(data$y2)+1,
            currentu=rep(0,nrow(data$y1)),gammatune=rep(0.000001,ncol(Za)),
            propa=1,propb=0.5,propx2=1/0.05,vx2=rep(10,nrow(Za)),
            x1propa=x1propa,x1propb=x1propb,
            currenttheta=rep(1,1),betaxtune=c(1,rep(0.01,ncol(Za)-1),1), 
            propax2=1,propbx2=0.5,
            currentlambda=0.5,propl1=5,propl2=23,
            currentdelta=1,propd1=6.6,propd2=3.6,alphatune=rep(0.0000001,ncol(data$Zb)),
            currentb=matrix(0,nrow=nrow(data$y1),ncol=2),btune=c(0.001,0.001),
            currentSigmab=diag(2)*0.00000001,currentzeta=sample(1:ncomp,nrow(data$y1),TRUE),
            currentpi=rep(1/ncomp,ncomp),currentm=apply(y1,1,max)+2,
            currentsigma2b=0.01,currentb2=matrix(0,nrow=nrow(data$y1),ncol=3),
            btune2=rep(0.001,3),currentSigmab2=diag(3)*0.01,currentb3=rep(0,nrow(data$y1)),
            currentphi=.89)

prior = list(mu0y2=rep(0,ncol(data$Za)+1),mu0x1=rep(0,ncol(Za)),mu0x2=rep(0,ncol(Za)+1),
             mu0a=rep(0,ncol(data$Zb)),V0y2=100*diag(ncol(data$Za)+1),V0x1=100*diag(ncol(Za)),
             V0x2=100*diag(ncol(Za)+1),V0a=100*diag(ncol(data$Zb)),a0eta=1,b0eta=1,
             a0x=1,b0x=1,a0y=1,b0y=1,
             a0theta=1,b0theta=1,
             a0l=1,b0l=1,
             a0delta=1,b0delta=1, d0=4, D0=diag(2),adirich=rep(1,ncomp),
             D02=diag(3))

mcmc1 = mcmc_2part_nci1(data=data,init=init,priors=prior,nrep=25000,burn=10000)
#mcmc2 = mcmc_2part_nci2(data=data,init=init,priors=prior,nrep=6000,burn=2000)
 mcmc3 = mcmc_2part_nci3(data=data,init=init,priors=prior,nrep=6000,burn=2000)
# mcmc4 = mcmc_2part_nci4(data=data,init=init,priors=prior,nrep=6000,burn=2000)

mcmc2b = mcmc_2part_nci2b(data=data,init=init,priors=prior,nrep=25000,burn=10000)
mcmc2c = mcmc_2part_nci2c(data=data,init=init,priors=prior,nrep=6000,burn=2000)
mcmc2d = mcmc_2part_nci2d(data=data,init=init,priors=prior,nrep=6000,burn=2000)

mcmc2c2 = mcmc_2part_nci2c2(data=data,init=init,priors=prior,nrep=6000,burn=2000)
mcmc2c3 = mcmc_2part_nci2c3(data=data,init=init,priors=prior,nrep=6000,burn=2000)

mcmc5 = mcmc_2part_nci5(data=data,init=init,priors=prior,nrep=6000,burn=2000)
mcmc5b = mcmc_2part_nci5b(data=data,init=init,priors=prior,nrep=20000,burn=10000)

mcmc6 = mcmc_2part_nci6(data=data,init=init,priors=prior,nrep=6000,burn=2000)

mcmc7 = mcmc_2part_nci7(data=data,init=init,priors=prior,nrep=100000,burn=20000)
mcmc7b = mcmc_2part_nci7b(data=data,init=init,priors=prior,nrep=20000,burn=10000)
mcmc7c = mcmc_2part_nci7c(data=data,init=init,priors=prior,nrep=20000,burn=10000)

mcmc8 = mcmc_2part_nci8(data=data,init=init,priors=prior,nrep=100000,burn=50000)
mcmc8b = mcmc_2part_nci8b(data=data,init=init,priors=prior,nrep=100000,burn=50000)

mcmc9 = mcmc_2part_nci9(data=data,init=init,priors=prior,nrep=6000,burn=2000)

#acceptance rates
apply(mcmc$gamma,2,function(x){return(length(unique(x))/length(x))}) 
apply(mcmc$betay,2,function(x){return(length(unique(x))/length(x))}) 
#apply(mcmc$alpha,2,function(x){return(length(unique(x))/length(x))}) 
plot(apply(mcmc$latentx1,2,function(x){return((length(unique(x))-1)/length(x))}))
#plot(apply(mcmc$latentx2,2,function(x){return((length(unique(x))-1)/length(x))}))
plot(apply(mcmc$b1,2,function(x){return((length(unique(x))-1)/length(x))}))
(length(unique(mcmc$eta))-1)/length(mcmc$eta)
#(length(unique(mcmc$theta))-1)/length(mcmc$theta)
(length(unique(mcmc$lambda))-1)/length(mcmc$lambda)
(length(unique(mcmc$delta))-1)/length(mcmc$delta)
(length(unique(mcmc$phi))-1)/length(mcmc$phi)

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
plot(mcmc$phi,type="l")


plot(mcmc$betay[,1],type="l")
plot(mcmc$betay[,2],type="l")
plot(mcmc$betay[,3],type="l")
plot(mcmc$betay[,4],type="l")
plot(mcmc$betay[,5],type="l")
plot(mcmc$betay[,6],type="l")
plot(mcmc$betay[,7],type="l")
plot(mcmc$betay[,8],type="l")
plot(mcmc$betay[,9],type="l")


plot(mcmc$sigma2b[,1],type="l")
plot(mcmc$sigma2b[,2],type="l")
plot(mcmc$sigma2b[,3],type="l")
plot(mcmc$corrb[,1],type="l")
plot(mcmc$corrb[,2],type="l")
plot(mcmc$corrb[,3],type="l")

plot(mcmc$sigma2y,type="l")



postcov <- do.call(cbind,mcmc)
postcov <- cor(with(mcmc,cbind(betay,gamma,phi,lambda,sigma2b)))
postcov[lower.tri(postcov)] <- 0
diag(postcov) <- 0
df <- expand.grid(x=1:ncol(postcov),y=1:ncol(postcov))
df$covar <- c(postcov)

ggplot(data=df,aes(x=as.factor(x),y=as.factor(y))) + geom_tile(aes(fill=covar)) + 
  #geom_text(aes(label=count)) + scale_fill_gradient(low = "white", high = "red") 
  scale_fill_gradient2(high="blue",low="red",mid="white")

nm <- data.frame(x=1:ncol(postcov),y=c(paste0("betay",1:ncol(mcmc$betay)),paste0("gamma",1:ncol(mcmc$gamma)),"phi","lambda",paste0("sigma2b",1:ncol(mcmc$sigma2b))))

coeffs <- cbind(mcmc$betay,mcmc$gamma)
cisb <- data.frame(t(apply(mcmc$betay[,1:8],2,quantile,probs=c(0.025,0.5,0.975))))
cisg <- data.frame(t(apply(mcmc$gamma,2,quantile,probs=c(0.025,0.5,0.975))))
cis <- rbind(cisb,cisg)
names(cis) <- c("q025","q50","q975")
cis$model <- c(rep("Weibull",nrow(cisb)),rep("GenPoisson",nrow(cisg)))
cis$name <- rep(c("Intercept","Age","Male","BMI","Smoke","Education","Black","Hispanic"),2)
cis$Significant <- "No"
cis$Significant[cis$q025 > 0 | cis$q975 < 0] <- "Yes"

ggplot(subset(cis,q50<4), aes(x=name,y=q50, ymin=q025, ymax=q975,colour=Significant)) + facet_grid(~model)+ geom_pointrange(fatten=5) + 
  geom_hline(yintercept=0) + coord_flip() + labs( y = 'Slopes') + theme_bw()

# mby <- melt(mcmc$betay[20000:100000,])
# mby$X2 <- factor(mby$X2,labels=c("Intercept","Age","Male","BMI","Smoke","Education","Black","Hispanic","log(Y1)"))
# ggplot(data=mby,aes(x=value))+geom_histogram() + facet_wrap(~X2,scales="free") + xlab(expression(beta)) + theme_bw()
# 
# may <- melt(mcmc$gamma[20000:100000,])
# may$X2 <- factor(may$X2,labels=c("Intercept","Age","Male","BMI","Smoke","Education","Black","Hispanic"))
# ggplot(data=may,aes(x=value))+geom_histogram() + facet_wrap(~X2,scales="free") + xlab(expression(gamma))+theme_bw()
#-------------------------------------------------------------

assessln <- pp_assess(mcmc1,data$Zb,1000,1,burn=10000)
# assessln2 <- pp_assess(mcmc2,data$Zb,200,2,burn=2000)
 assessln3 <- pp_assess(mcmc3,data$Zb,400,3,y1,y2,burn=2000)
# assessln4 <- pp_assess(mcmc4,data$Zb,200,4)
assessln2b <- pp_assess(mcmc2b,data$Zb,1000,"2b",burn=10000)
assessln2c <- pp_assess(mcmc2c,data$Zb,400,"2c",y1,y2,burn=2000)
assessln2d <- pp_assess(mcmc2d,data$Zb,400,"2d",y1,y2,burn=2000)

assessln2c2 <- pp_assess(mcmc2c2,data$Zb,400,"2c2",y1,y2,burn=2000)
assessln2c3 <- pp_assess(mcmc2c3,data$Zb,400,"2c3",y1,y2,burn=2000)

assessln5 <- pp_assess(mcmc5,data$Zb,400,"5",y1,y2,burn=2000)
assessln5b <- pp_assess(mcmc5b,data$Zb,400,"5b",y1,y2,burn=10000)

assessln6 <- pp_assess(mcmc6,data$Zb,400,"6",y1,y2,burn=10000)

assessln7 <- pp_assess(mcmc7,data$Zb,1000,"7",y1,y2,weights,burn=20000)
assessln7b <- pp_assess(mcmc7b,data$Zb,400,"7b",y1,y2,weights,burn=10000)
assessln8 <- pp_assess(mcmc8,data$Zb,2000,"8",y1,y2,burn=50000)
assessln8b <- pp_assess(mcmc8b,data$Zb,2000,"8b",y1,y2,burn=50000)

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
y2median <- median(c(y2[y2>0]))
y2mean <- mean(c(y2[y2>0]))
y2var <- var(c(y2[y2>0]))
y2q15 <- quantile(c(y2[y2>0]),probs=c(0.20))
y2q25 <- quantile(c(y2[y2>0]),probs=c(0.25))
y2q35 <- quantile(c(y2[y2>0]),probs=c(0.35))
y2q90 <- quantile(c(y2[y2>0]),probs=c(0.9))
y2q95 <- quantile(c(y2[y2>0]),probs=c(0.95))
y2q30 <- quantile(c(y2[y2>0]),probs=c(0.3))

y100 <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
y110 <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
y101 <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
y120 <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
y102 <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
y112 <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
y121 <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
y111 <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
y122 <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))

probs <- c(y100,y110,y120,y101,y102,y111,y112,y121,y122)/sum(c(y100,y110,y120,y101,y102,y111,y112,y121,y122))


y2daydiff <- mean(y2[,1]-y2[,2])

assess=assessln$out
obs <- c(median(assess$y100),median(assess$y110),median(assess$y120),median(assess$y101),
         median(assess$y102),median(assess$y111),median(assess$y112),median(assess$y121),
         median(assess$y122))
chisq.test(obs,p=probs)
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
q9b2 <- qplot(x=assess$y2q30) + geom_vline(xintercept=y2q30,colour="red") + theme_bw()
q9c <- qplot(x=assess$y2q35) + geom_vline(xintercept=y2q35,colour="red") + theme_bw()
q10 <- qplot(x=assess$y2q90) + geom_vline(xintercept=y2q90,colour="red") + theme_bw()
q11 <- qplot(x=assess$y2q95) + geom_vline(xintercept=y2q95,colour="red") + theme_bw()
q12 <- qplot(x=assess$y2mean) + geom_vline(xintercept=y2mean,colour="red") + theme_bw()
q13 <- qplot(x=assess$y2var) + geom_vline(xintercept=y2var,colour="red") + theme_bw()

grid.arrange(q1,q2,q2b,q2c,q3,q4,q6,q7b,q8,q9a,q9b,q9c,q10,q11,q12,q13)

comply <- data.frame(No=assess$pcomply,Yes=assess$pcomply2)
mcomply <- melt(comply)
names(mcomply)[1] <- "Weight"
ci <- quantile(comply$No,probs=c(0.025,0.975))
wci <- quantile(comply$Yes,probs=c(0.025,0.975))
qplot(data=mcomply,x=value,fill=Weight,geom="density",alpha=I(0.5)) + geom_vline(xintercept=ci,linetype=2,colour="red") + 
  geom_vline(xintercept=wci, linetype=2, colour="blue") + ggtitle("Distribution of Compliance Rates to PAG for Entire Population") + theme_bw()
