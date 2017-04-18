library(MASS)
library(MCMCpack)
library(mvtnorm)

pams <- read.csv("C:\\Users\\dcries\\github\\pams\\FinalPAMSDataSetNew.csv")
pams <- pams[-which(is.na(pams$ModPAR)|is.na(pams$Age)|is.na(pams$BMI)|is.na(pams$Smoker)),]


log_q <- function(y,b,p,mu,sigma2,Sigmab){
  ll1 <- dmvnorm(b,c(0,0),Sigmab,log=TRUE)
  if(y==0){
    ll2 <- log(1-p)
  }
  else{
    ll2 <- log(p) + dlnorm(y,mu,sqrt(sigma2),log=TRUE)
  }
  return(ll1+ll2)
}

mcmc_ln <- function(data,init,prior,nreps,burn=1000){
  #data
  Xa    <- data$Xa
  Xb <- data$Xb
  y <- data$y
  n <- nrow(Xa)
  z <- rep(1,n)
  z[y==0] <- 0
  ind <- y>0
  yp <- y[ind]
  Xbp <- Xb[ind,]
  
  #initial values
  currentbeta <- init$currentbeta
  currentalpha <- init$currentalpha
  currentsigma2 <- init$currentsigma2
  currentSigmab <- init$currentSigmab
  currentb <- init$currentb
  currentu <- rep(0,n)
  
  #priors
  mu0a <- prior$mu0a
  mu0b <- prior$mu0b
  V0a <- prior$V0a
  V0b <- prior$V0b
  a0 <- prior$a0
  b0 <- prior$b0
  nu0 <- prior$nu0
  D0 <- prior$D0
  
  #allocate storae
  beta <- matrix(0,ncol=ncol(Xb),nrow=nreps)
  alpha <- matrix(0,ncol=ncol(Xa),nrow=nreps)
  sigma2 <- rep(0,nreps)
  #u <- matrix(0,)
  sigma2b1 <- rep(0,nreps)
  sigma2b2 <- rep(0,nreps)
  corrb <- rep(0,nreps)
  b1 <- matrix(0,ncol=n,nrow=nreps)
  b2 <- matrix(0,ncol=n,nrow=nreps)
  
  tune <- array(0.01*diag(2),dim=c(2,2,n))
  
  for(i in 1:nreps){
    for(j in 1:n){
      rn <- rnorm(1,Xa[j,]%*%currentalpha+currentb[j,1],1)
      currentu[j] <- ifelse(y[j]==0,-abs(rn),abs(rn))
    }
    
    #sample alphas
    Va <- solve(solve(V0a)+t(Xa)%*%Xa)
    mua <- Va%*%(solve(V0a)%*%mu0a + t(Xa)%*%(currentu-currentb[,1]))
    currentalpha <- mvrnorm(1,mua,Va)
    currentp <- pnorm(Xa%*%currentalpha+currentb[,1])
    
    #sample betas
    Vb <- solve(solve(V0b)+t(Xbp)%*%Xbp/currentsigma2)
    mub <- Vb%*%(solve(V0b)%*%mu0b + t(Xbp)%*%(log(yp)-currentb[ind,2])/currentsigma2)
    currentbeta <- mvrnorm(1,mub,Vb)
    currentmu <- exp(Xb%*%currentbeta+currentb[,2])
    
    #sample sigma2
   currentsigma2 <- rinvgamma(1,n/2+a0,(b0+0.5*sum((log(yp)-Xbp%*%currentbeta-currentb[ind,2])^2)))

    
    #sample covaraicne matrix for re
    currentSigmab <- riwish(nu0+n,D0+t(currentb)%*%currentb)
    
    for(j in 1:n){
      propb <- mvrnorm(1,currentb[j,],2.88*tune[,,j])
      lacceptprob <- log_q(y[j],propb,currentp[j],currentmu[j],currentsigma2,currentSigmab) - 
        log_q(y[j],currentb[j,],currentp[j],currentmu[j],currentsigma2,currentSigmab)
      
      if(lacceptprob > log(runif(1))){
        currentb[j,] <- propb
      }
      if((i>20)&(i<burn)&(i%%20==0)){
        tune[,,j] <- cov(cbind(b1[1:(i-1),j],b2[1:(i-1),j]))
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
    
    if(i%%100==0){
      print(i)
    }
  }
  params <- data.frame(cbind(alpha,beta,sigma2,sigma2b1,sigma2b2,corrb))
  names(params) <- c(paste0("alpha_",0:(ncol(Xa)-1)),
                paste0("beta_",0:(ncol(Xb)-1)),
                "sigma2","sigma2b1","sigma2b2","corrb")
  
  out <- list(params=params,b1=b1,b2=b2)
  return(out)
} 

#pams <- pams[1:500,]
data <- list(Xa=with(pams,model.matrix(~Age+BMI+Gender+Smoker)),
             Xb=with(pams,model.matrix(~Age+BMI+Gender+Smoker)),
                      y=pams$ModPAR)
init <- list(currentbeta=rep(0,ncol(data$Xb)),
                currentalpha=rep(0,ncol(data$Xa)),
                currentsigma2=1,
                currentSigmab=diag(2),
                currentb=matrix(0,ncol=2,nrow=nrow(pams)))
prior <- list(mu0a=rep(0,ncol(data$Xa)),
              mu0b=rep(0,ncol(data$Xb)),
              V0a=10*diag(ncol(data$Xa)),
              V0b=10*diag(ncol(data$Xb)),
              a0=1,
              b0=1,
              nu0=3,
              D0=diag(2))
                      
out <- mcmc_ln(data,init,prior,10000)
