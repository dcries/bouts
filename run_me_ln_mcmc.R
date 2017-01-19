Rcpp::sourceCpp('C:/Users/dcries/github/bouts/mcmc_me_ln.cpp')
Rcpp::sourceCpp('/home/danny/Documents/github/bouts/mcmc_me_ln.cpp')
Rcpp::sourceCpp('/home/danny/Documents/github/bouts/ppred_check.cpp')


pams <- read.csv("C:\\Users\\dcries\\github\\pams\\FinalPAMSDataSetNew.csv")
pams <- read.csv("/home/danny/Documents/github/pams/FinalPAMSDataSetNew.csv")

pams <- pams[-which(is.na(pams$ModPAR)|is.na(pams$Age)|is.na(pams$BMI)|is.na(pams$Smoker)),]
p1 <- pams[pams$Trial==1,]
p2 <- pams[pams$Trial==2,]

pc <- merge(p1,p2,by="id")
pcy <- cbind(pc$ModPAR.x,pc$ModPAR.y)
pcw <- cbind(pc$Moderatev52.x,pc$Moderatev52.y)

pcy[pcy<10] <- 0
pcw[pcw<10] <- 0


data <- list(Za=with(pc,model.matrix(~Age.y+BMI.y+Gender.y+Smoker.y)),
             Zb=with(pc,model.matrix(~Age.y+BMI.y+Gender.y+Smoker.y)),
             y=pcy,w=pcw)

init <- list(currentbeta=rep(0,ncol(data$Zb)),
             currentalpha=rep(0,ncol(data$Za)),
             currentsigma2=1,
             currentSigmab=diag(2),
             currentb=matrix(0,ncol=2,nrow=nrow(pc)),
             currentsigma2w=1,
             currentsigma2y=1,
             currentgamma=1,
             currentx=rowMeans(pcw),
             tune=rep(0.001,2))

prior <- list(mu0a=rep(0,ncol(data$Za)),
              mu0b=rep(0,ncol(data$Zb)),
              mu0g=0,
              V0a=100*diag(ncol(data$Za)),
              V0b=100*diag(ncol(data$Zb)),
              V0g=100,
              a0=1,b0=1,a0w=1,b0w=1,a0y=1,b0y=1,
              nu0=3,
              D0=diag(2))

out <- mcmc_me_ln(data,init,prior,10000,burn=2000,p0=0.85)


check <- ppred_ln(out,data,1000)
summary(check$zeros); apply(pcw,2,function(x){sum(x==0)})
summary(check$maxw); apply(pcw,2,max)
summary(check$maxy); apply(pcy,2,max)
summary(check$complyw); apply(pcw,2,function(x){sum(x>450/7)/length(x)})
summary(check$complyy); apply(pcy,2,function(x){sum(x>450/7)/length(x)})


out$accept0/out$prop0;out$accept1/out$prop1
plot(apply(out$b1,2,function(x){(length(unique(x))-1)/length(x)}))
prob <- pnorm(data$Za%*%colMeans(out$alpha)+colMeans(out$b1))
plot(prob)
#prob vs modpar
plot((sqrt(out$sigma2)))
plot((sqrt(out$sigma2w)))
plot((sqrt(out$sigma2y)))
plot((sqrt(out$sigma2b1)))
plot((sqrt(out$sigma2b2)))
plot(((out$corrb)))


plot(pcw[,1],prob)
hist(out$latentx[,which.max(apply(out$latentx,2,function(x){sum(x==0)}))])
hist(out$latentx[,which.max(prob)])


mu <- data$Zb%*%colMeans(out$beta) + colMeans(out$b2)
hist(exp(mu)) #hmmmm

#residuals for lognormal model, obs > 0
residw1 <- mu[pcw[,1]>0] - log(pcw[pcw[,1]>0,1])
residw2 <- mu[pcw[,2]>0] - log(pcw[pcw[,2]>0,2])
residy1 <- mu[pcy[,1]>0] - log(pcy[pcy[,1]>0,1])
residy2 <- mu[pcy[,2]>0] - log(pcy[pcy[,2]>0,2])

#residuals for latentx obs >0
truex <- (colMeans(out$latentx))
residx <- mu - log(truex)

residw <- c(residw1,residw2)
residy <- c(residy1,residy2)

truew <- c(pcw[pcw[,1]>0,1],pcw[pcw[,2]>0,2])
truey <- c(mu[pcy[,1]>0],mu[pcy[,2]>0])


plot(log(truew),residw)
plot(truey,residy)
plot(truex,residx)

#qq plots
qqnorm(residw);qqline(residw);shapiro.test(residw)
qqnorm(residy);qqline(residy);shapiro.test(residy)
qqnorm(residx);qqline(residx);shapiro.test(residx)

#not sure what these would even tell us
qqnorm(colMeans(out$b1))
qqnorm(colMeans(out$b2))
