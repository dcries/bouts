library(MASS)
library(MCMCpack)
library(mvtnorm)

pams <- read.csv("C:\\Users\\dcries\\github\\pams\\FinalPAMSDataSetNew.csv")
pams <- pams[-which(is.na(pams$ModPAR)|is.na(pams$Age)|is.na(pams$BMI)|is.na(pams$Smoker)),]


log_qu <- function(y,b,p,mu,sigma2,Sigmab){
  ll1 <- dmvnorm(b,c(0,0),Sigmab,log=TRUE)
  if(y==0){
    ll2 <- log(1-p)
  }
  else{
    ll2 <- log(p) + dlnorm(y,mu,sqrt(sigma2),log=TRUE)
  }
  #ll2 <- dlnorm(y,mu,sqrt(sigma2),log=TRUE)
  return(ll1+ll2)
}

log_qx <- function(x,w,y,p,mu,gamma,sigma2,sigma2w,sigma2y){ 
  if(x==0){
    ll1 <- log(1-p)
  }
  else{
    ll1 <- log(p) + dlnorm(x,mu,sqrt(sigma2),log=TRUE)
  }
  #ll1 <- dlnorm(x,mu,sqrt(sigma2),log=TRUE)
  ll3 <- sum(dnorm(w,x,sqrt(sigma2w),log=TRUE) + dnorm(y,x*gamma,sqrt(sigma2y),log=TRUE))
  
  return(ll1+ll3)
}  

log_gx <- function(x,p0,scl){
  if(x==0){
    ll <- log(1-p0)
  }
  else{
    ll <- log(p0) + dexp(x,0.0095,log=TRUE)#2*dt(x,ncp=10000,df=20,log=TRUE)
  }
  #ll <- dexp(x,0.0095,log=TRUE)#2*dt(x,ncp=10000,df=20,log=TRUE)
  return(ll)
}

mcmc_ln <- function(data,init,prior,nreps,burn=1000){
  #!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!
  #non adaptive right now
  p0 <- 0.85
  
  #data
  Za    <- data$Za
  Zb <- data$Zb
  y <- data$y
  w <- data$w
  n <- nrow(Za)
  
  nr <- ncol(y)
  ybar <- rowMeans(y)
  scl <- 2#quantile(rowMeans(w),probs=0.75)
  
  #initial values
  currentbeta <- init$currentbeta
  currentalpha <- init$currentalpha
  currentsigma2 <- init$currentsigma2
  currentSigmab <- init$currentSigmab
  currentb <- init$currentb
  currentsigma2w <- init$currentsigma2w
  currentsigma2y <- init$currentsigma2y
  currentgamma <- init$currentgamma
  currentx <- init$currentx
  currentu <- rep(0,n)
  
  #priors
  mu0a <- prior$mu0a
  mu0b <- prior$mu0b
  mu0g <- prior$mu0g
  V0a <- prior$V0a
  V0b <- prior$V0b
  V0g <- prior$V0g
  a0 <- prior$a0
  b0 <- prior$b0
  a0w <- prior$a0w
  b0w <- prior$b0w
  a0y <- prior$a0y
  b0y <- prior$b0y
  nu0 <- prior$nu0
  D0 <- prior$D0
  
  #allocate storae
  beta <- matrix(0,ncol=ncol(Zb),nrow=nreps)
  alpha <- matrix(0,ncol=ncol(Za),nrow=nreps)
  sigma2 <- rep(0,nreps)
  #u <- matrix(0,)
  sigma2b1 <- rep(0,nreps)
  sigma2b2 <- rep(0,nreps)
  corrb <- rep(0,nreps)
  b1 <- matrix(0,ncol=n,nrow=nreps)
  b2 <- matrix(0,ncol=n,nrow=nreps)
  gamma <- rep(0,nreps)
  sigma2w <- rep(0,nreps)
  sigma2y <- rep(0,nreps)
  latentx <- matrix(0,ncol=n,nrow=nreps)
  
  
  tune <- array(0.001*diag(2),dim=c(2,2,n))
  
  #to check acceptance rates
  prop0 <- 0
  accept0 <- 0
  prop1 <- 0
  accept1 <- 0
  
  for(i in 1:nreps){
    
    ind <- currentx>0
    currentxp <- currentx[ind]
    Zbp <- Zb[ind,]
    
    for(j in 1:n){
      rn <- rnorm(1,Za[j,]%*%currentalpha+currentb[j,1],1)
      currentu[j] <- ifelse(currentx[j]==0,-abs(rn),abs(rn))
    }
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #rn <- rnorm(1,Za%*%currentalpha+currentb[,1],1)
    #currentu <- ifelse(currentx==0,-abs(rn),abs(rn))
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    #sample alphas
    Va <- solve(solve(V0a)+t(Za)%*%Za)
    mua <- Va%*%(solve(V0a)%*%mu0a + t(Za)%*%(currentu-currentb[,1]))
    currentalpha <- mvrnorm(1,mua,Va)
    currentp <- pnorm(Za%*%currentalpha+currentb[,1])
    
    #sample betas
    Vb <- solve(solve(V0b)+t(Zbp)%*%Zbp/currentsigma2)
    mub <- Vb%*%(solve(V0b)%*%mu0b + t(Zbp)%*%(log(currentxp)-currentb[ind,2])/currentsigma2)
    currentbeta <- mvrnorm(1,mub,Vb)
    currentmu <- exp(Zb%*%currentbeta+currentb[,2])
    
    #sample sigma2
    currentsigma2 <- rinvgamma(1,n/2+a0,(b0+0.5*sum((log(currentxp)-Zbp%*%currentbeta-currentb[ind,2])^2)))
    
    
    #sample covaraicne matrix for re
    currentSigmab <- riwish(nu0+n,D0+t(currentb)%*%currentb)
    
    for(j in 1:n){
      propb <- mvrnorm(1,currentb[j,],2.88*tune[,,j])
      propp <- pnorm(Za[j,]%*%currentalpha+propb[1])
      propmu <- exp(Zb[j,]%*%currentbeta+propb[2])
      
      lacceptprobu <- log_qu(currentx[j],propb,propp,log(propmu),currentsigma2,currentSigmab) - 
        log_qu(currentx[j],currentb[j,],currentp[j],log(currentmu[j]),currentsigma2,currentSigmab)
      
      if(lacceptprobu > log(runif(1))){
        currentb[j,] <- propb
      }
      if((i>20)&(i<burn)&(i%%20==0)){
        tune[,,j] <- cov(cbind(b1[1:(i-1),j],b2[1:(i-1),j]))
      }
      
    } 
    
    #sample sigma2y sigma2w
    wsum <- 0
    ysum <- 0
    for(j in 1:nr){
      wsum <- wsum + sum((w[,j]-currentx)^2)
      ysum <- ysum + sum((y[,j]-currentgamma*currentx)^2)
    }
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #wsum <- sum((w-currentx)^2)
    #ysum <- sum((y-currentgamma*currentx)^2)
    
    currentsigma2w <- rinvgamma(1,nr*n/2+a0w,b0w+0.5*wsum)
    currentsigma2y <- rinvgamma(1,nr*n/2+a0y,b0y+0.5*ysum)
    
    #sample gamma
    Vg <- solve(t(currentx)%*%currentx/currentsigma2y + 1/V0g)
    mg <- Vg*(mu0g/V0g+t(currentx)%*%ybar/currentsigma2y)
    currentgamma <- rnorm(1,mg,Vg)
    
    for(j in 1:n){
      r <- runif(1)
      #propx <- ifelse(r<(1-p0),0,rexp(1,.0095))#abs(rt(1,ncp=10,df=20)))
      if(r<(1-p0)){
        propx <- 0
        prop0 <- prop0+1
      }
      else{
        propx <- rexp(1,.0095)
        prop1 <- prop1+1
      }
      lacceptprobx <- log_qx(propx,w[j,],y[j,],currentp[j],log(currentmu[j]),currentgamma,currentsigma2,currentsigma2w,currentsigma2y) + log_gx(currentx[j],p0,scl) -
        log_qx(currentx[j],w[j,],y[j,],currentp[j],log(currentmu[j]),currentgamma,currentsigma2,currentsigma2w,currentsigma2y) - log_gx(propx,p0,scl)
      
      if(lacceptprobx > log(runif(1))){
        currentx[j] <- propx
        if(propx==0){
          accept0 <- accept0+1
        } 
        else{
          accept1 <- accept1+1
        }
      }
    }
    
    
    
    beta[i,] <- currentbeta
    alpha[i,] <- currentalpha
    sigma2[i] <- currentsigma2
    #u <- matrix(0,)
    sigma2b1[i] <- currentSigmab[1,1]
    sigma2b2[i] <- currentSigmab[2,2]
    corrb[i] <- currentSigmab[1,2]/sqrt(currentSigmab[1,1]*currentSigmab[2,2])
    b1[i,] <- t(currentb[,1])
    b2[i,] <- t(currentb[,2])
    sigma2w[i] <- currentsigma2w
    sigma2y[i] <- currentsigma2y
    gamma[i] <- currentgamma
    latentx[i,] <- currentx
    
    if(i%%100==0){
      print(i)
    }
  }
  params <- data.frame(cbind(alpha,beta,gamma,sigma2,sigma2w,sigma2y,
                             sigma2b1,sigma2b2,corrb))
  names(params) <- c(paste0("alpha_",0:(ncol(Za)-1)),
                     paste0("beta_",0:(ncol(Zb)-1)),
                     "gamma","sigma2","sigma2w","sigma2y",
                     "sigma2b1","sigma2b2","corrb")
  
  out <- list(params=params,b1=b1,b2=b2,latentx=latentx)
  return(out)
} 

#pams <- pams[1:500,]
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
             currentx=rowMeans(pcw))

prior <- list(mu0a=rep(0,ncol(data$Za)),
              mu0b=rep(0,ncol(data$Zb)),
              mu0g=0,
              V0a=10*diag(ncol(data$Za)),
              V0b=10*diag(ncol(data$Zb)),
              V0g=10,
              a0=1,b0=1,a0w=0,b0w=0,a0y=0,b0y=0,
              nu0=3,
              D0=diag(2))

out <- mcmc_ln(data,init,prior,10000,burn=1000)
